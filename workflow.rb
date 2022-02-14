require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/MODULE'

Workflow.require_workflow "HTS"
Workflow.require_workflow "Sequence"
Workflow.require_workflow "PCAWG"
module PCAWGPilot48
  extend Workflow

  SAMPLES = Rbbt.data.donors.list - %w(DO50311)

  input :vcf, :path, "VCF file"
  input :bed, :path, "BED file"
  extension :vcf
  task :subset_by_bed do |vcf,bed|
    output = file('output/donor')
    Open.mkdir file('output')
    #CMD.cmd("vcftools --vcf #{vcf} --bed #{bed}  --out #{output} --recode --keep-INFO-all")
    #Open.mv file('output').glob("*.vcf").first, self.tmp_path
    if Misc.is_filename?(vcf)
      CMD.cmd("bcftools view --targets-file  #{bed} #{vcf}  -o #{self.tmp_path}")
    else
      TmpFile.with_file(vcf) do |tvcf|
        CMD.cmd("bcftools view --targets-file  #{bed} #{tvcf}  -o #{self.tmp_path}")
      end
    end
    nil
  end

  dep :subset_by_bed
  input :mutation_type, :select, "Mutations to consider", :SNV, :select_options => %w(SNV indel both)
  extension :vcf
  task :subset_by_type => :text do |mutation_type|
    TSV.traverse step(:subset_by_bed), :type => :array, :into => :stream do |line|
      next line if line =~ /^#/
      parts = line.split("\t")

      ref, alt = parts.values_at 3, 4

      snv = true
      snv = false unless %w(A C T G).include? ref
      snv = false unless %w(A C T G).include? alt

      case mutation_type.to_s.downcase
      when 'snv'
        next unless snv
      when 'indel'
        next if snv
      when 'both'
      else
        raise ParameterException, "Unknown mutation type: #{mutation_type}"
      end

      line
    end
  end

  dep :subset_by_type, :compute => :produce
  input :maf, :float, "Minimun maf", 0.05
  input :pass_only, :boolean, "Condider only PASS variants", true
  extension :vcf
  task :subset_by_af => :text do |maf,pass_only|
    vcf = step(:subset_by_type).join.path

    fields = TSV.parse_header(vcf).fields

    tumor_sample = begin
                     CMD.cmd("grep 'tumor_sample=' '#{vcf}'").read.strip.split("=").last
                   rescue
                     if TSV.parse_header(vcf).fields.include? "TUMOR"
                       Log.warn "Could not find tumor_sample field in input VCF, using TUMOR"
                       "TUMOR"
                     else
                       Log.warn "Could not find tumor_sample field in input VCF, using last field"
                       TSV.parse_header(vcf).fields.last
                     end
                   end

    tumor_pos = fields.index tumor_sample

    TSV.traverse vcf, :type => :array, :into => :stream do |line|
      next line if line =~ /^##/
      if line =~ /^#CHR/
        format_line = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' << "\n"
        next format_line + line 
      end
      parts = line.split("\t")

      next if pass_only  && (parts[6] != "PASS")

      format = parts[8].split(":")
      donor = parts[tumor_pos + 1].split(":")
      hash = Misc.zip2hash(format, donor)

      if hash["AF"]
        next if hash["AF"].to_f <= maf
      elsif hash["TIR"]
        tar = hash["TAR"].split(",").first.to_f
        tir = hash["TIR"].split(",").first.to_f
        af = tir.to_f / (tar + tir)
        next if af <= maf
      elsif hash["AU"]
        a,c,t,g = hash.values_at(*%w(AU CU TU GU)).collect{|v| v.to_i}
        sum = a+c+t+g
        var = case parts[4]
              when "A"
                a
              when "C"
                c
              when "T"
                t
              when "G"
                g
              else
                raise parts[4]
              end
        af = var.to_f / sum
        next if af <= maf
      elsif hash["AD"]
        parts = hash["AD"].split(",")
        ref, alt, *rest = parts
        total = parts.inject(0){|acc,e| acc += e.to_i}
        af = alt.to_f / total
        next if af <= maf
      elsif hash["PM"]
        hash["PM"].to_f
      elsif hash["PR"]
        total = Misc.sum(hash.values_at("PR", "NR").collect(&:to_f))
        support = Misc.sum(hash.values_at("PU", "NU").collect(&:to_f))
        support.to_f / total.to_f
      else
        raise "No entry with AF in VCF file"
      end

      if hash.include? "GT"
        if hash["GT"] == "./."
          parts[tumor_pos + 1] = parts[tumor_pos + 1].sub("./.", "1/0")
          parts * "\t"
        else
          line
        end
      else
        parts[8] += ":GT"
        parts[tumor_pos + 1] += ":1/0"
        parts * "\t"
      end
    end
  end

  input :bed_type, :select, "BED type", :capture, :select_options => %w(capture 50 100 150 200)
  extension :bed
  task :bed_file => :path do |bed_type|
    donor ||= clean_name
    file = if bed_type.to_s == 'capture'
             tsv = Rbbt.data["release_may2016.v1.4.tsv"].tsv(:header_hash => '', :fields => ["tumor_wgs_bwa_alignment_gnos_id"], :key_field => "icgc_donor_id", :type => :single)
             id = tsv[donor]
             num = Rbbt.data.BEDs.ICGC64_validation_beds["map.txt"].tsv(:key_field => 1, :fields => [0], :type => :single)[id]
             Rbbt.data.capture_beds["Array#{num}.bed"].find
           else
             Rbbt.data.computed_beds[bed_type].glob(donor + "-DS*").first 
           end
    Open.cp file, self.tmp_path
    nil
  end

  dep :bed_file
  input :ds_caller, :select, "Caller used for DS", :mutect2, :select_options => %w(mutect2 strelka muse)
  input :ds_system, :select, "System used for DS", :rbbt, :select_options => %w(rbbt rbbt.1 ARGO-aln ARGO)
  extension :vcf
  dep_task :ds_vcf, PCAWGPilot48, :subset_by_af, :bed => :bed_file, :vcf => :placeholder  do |donor,options|
    file = case options[:ds_system].to_s
           when 'rbbt'
             Rbbt.gold_standard[options[:ds_caller].to_s][donor + "-DS-orig.vcf"].find
           when 'rbbt.1'
             Rbbt.gold_standard[options[:ds_caller].to_s][donor + "-DS-rbbt.vcf"].find
           when 'ARGO-aln'
             Rbbt.gold_standard[options[:ds_caller].to_s][donor + "-DS.vcf"].find
           when 'ARGO'
             Rbbt.gold_standard[options[:ds_caller].to_s][donor + "-DS-ARGO.vcf"].find
           end
    {:inputs => options.merge(:vcf => file), :jobname => donor + "-DS"}
  end

  dep :bed_file
  input :wgs_caller, :select, "Caller used for WGS", :mutect2, :select_options => %w(mutect2 strelka muse)
  extension :vcf
  dep_task :wgs_vcf, PCAWGPilot48, :subset_by_af, :bed => :bed_file, :vcf => :placeholder  do |donor,options|
    file = Rbbt.result[options[:wgs_caller].to_s][donor + "-WGS.vcf"].find
    {:inputs => options.merge(:vcf => file), :jobname => donor + "-WGS"}
  end

  dep :wgs_vcf
  dep :ds_vcf
  input :donor, :select, "Donor to compare to", nil, :select_options => SAMPLES, :jobname => true
  dep_task :donor_vcfeval, HTS, :vcfeval, :input_vcf => :wgs_vcf, :truth_vcf => :ds_vcf, :reference => 'hg38'

  input :vcf_file, :file, "VCF file to evaluate", nil, :required => true
  input :donor, :select, "Donor to compare to", nil, :select_options => SAMPLES, :jobname => true
  dep :bed_file, :jobname => :donor
  dep :subset_by_af, :vcf => :vcf_file, :bed => :bed_file, :compute => :produce, :jobname => :donor
  dep_task :evaluate, PCAWGPilot48, :donor_vcfeval, :wgs_caller => :placeholder, "PCAWGPilot48#wgs_vcf" => :subset_by_af, :jobname => :donor

end

require 'tasks/genomic_mutations'
require 'tasks/figures'

#require 'rbbt/knowledge_base/MODULE'
#require 'rbbt/entity/MODULE'

