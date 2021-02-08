require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/MODULE'

Workflow.require_workflow "HTS"
Workflow.require_workflow "Sequence"
Workflow.require_workflow "PCAWG"
module PCAWGPilot
  extend Workflow

  SAMPLES = Rbbt.data.samples.list

  input :vcf, :file, "VCF file"
  input :bed, :file, "BED file"
  task :subset_by_bed do |vcf,bed|
    output = file('output/sample')
    Open.mkdir file('output')
    #CMD.cmd("vcftools --vcf #{vcf} --bed #{bed}  --out #{output} --recode --keep-INFO-all")
    #Open.mv file('output').glob("*.vcf").first, self.tmp_path
    CMD.cmd("bcftools view --targets-file  #{bed} #{vcf}  -o #{self.tmp_path}")
    nil
  end

  dep :subset_by_bed
  input :maf, :float, "Minimun maf", 0.05
  task :subset_by_af => :text do |maf|
    vcf = step(:subset_by_bed).path

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

      format = parts[8].split(":")
      sample = parts[tumor_pos + 1].split(":")
      hash = Misc.zip2hash(format, sample)
      if hash["AF"]
        next if hash["AF"].to_f <= maf
      elsif hash["TIR"]
        tar = hash["TAR"].split(",").first.to_f
        tir = hash["TIR"].split(",").first.to_f
        af = tir / (tar + tir)
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
      end

      if hash.include? "GT"
        line
      else
        parts[8] += ":GT"
        parts[tumor_pos + 1] += ":1/0"
        parts * "\t"
      end
    end
  end

  input :ds, :file
  input :wgs, :file
  dep :subset_by_af, :vcf => :ds, :jobname => "DS"
  dep :subset_by_af, :vcf => :wgs, :jobname => "WGS"
  dep_task :vcfeval, HTS, :vcfeval, :input_vcf => :placeholder, :truth_vcf => :placeholder, :reference => 'hg38' do |jobname,options, dependencies|
    ds, wgs = dependencies.flatten.values_at 0, 1
    options[:input_vcf] = wgs
    options[:truth_vcf] = ds
    {:inputs => options, :jobname => jobname}
  end

  dep :vcfeval
  dep_task :fp_positions, Sequence, :genomic_mutations, :vcf_file => :placeholder do |jobname,options,dependencies|
    vcfeval = dependencies.flatten.first
    options = options.merge(:vcf_file => vcfeval.file('output/fp.vcf.gz'))
    {:inputs => options, :jobname => jobname}
  end

  dep :vcfeval
  dep_task :fn_positions, Sequence, :genomic_mutations, :vcf_file => :placeholder do |jobname,options,dependencies|
    vcfeval = dependencies.flatten.first
    options = options.merge(:vcf_file => vcfeval.file('output/fn.vcf.gz'))
    {:inputs => options, :jobname => jobname}
  end

  input :bed_type, :select, "BED type", :capture, :select_options => %w(capture 50 100 150 200)
  input :ds_caller, :select, "Caller used for DS", :mutect2, :select_options => %w(mutect2 strelka)
  input :wgs_caller, :select, "Caller used for WGS", :mutect2, :select_options => %w(mutect2 strelka)
  dep_task :sample_vcfeval, PCAWGPilot, :vcfeval, :ds => :placeholder, :wgs => :placeholder do |sample,options,dependencies|
    ds_caller = options[:ds_caller] || 'mutect2'
    wgs_caller = options[:wgs_caller] || 'mutect2'
    options[:ds] = Rbbt.data.vcf[ds_caller.to_s][sample + "-DS.vcf"].find
    options[:wgs] = Rbbt.data.vcf[wgs_caller][sample + "-WGS.vcf"].find

    if options[:bed_type].to_s == 'capture'
      tsv = Rbbt.data["release_may2016.v1.4.tsv"].tsv(:header_hash => '', :fields => ["tumor_wgs_bwa_alignment_gnos_id"], :key_field => "icgc_donor_id", :type => :single)
      id = tsv[sample]
      num = Rbbt.data.BEDs.ICGC64_validation_beds["map.txt"].tsv(:key_field => 1, :fields => [0], :type => :single)[id]
      options[:bed] ||= Rbbt.data.capture_beds["Array#{num}.bed"].find
    else
      options[:bed] ||= Rbbt.data.computed_beds[options[:bed_type]].glob(sample + "-DS*").first 
    end

    {:inputs => options}
  end

  dep :sample_vcfeval, :compute => [:bootstrap, 20], :maf => :placeholder, :bed_type => :placeholder do |jobname,options,dependencies|
    Rbbt.data.samples.list.collect do |sample|
      [0.1, 0.05, 0.01, 0.005, 0.001].collect do |maf|
        %w(capture 50 100 150 200).collect do |bed_type|
          {:inputs => options.merge(:maf => maf, :bed_type => bed_type), :jobname => sample}
        end
      end
    end.flatten
  end
  task :explore_values => :tsv do
    all = nil
    dependencies.each do |dep|
      tsv = dep.load
      maf = dep.recursive_inputs[:maf]
      bed_type = dep.recursive_inputs[:bed_type]
      sample = dep.clean_name

      values = tsv[tsv.keys.first]

      tsv.keys.each{|k| tsv.delete k}
      key = [sample, bed_type, maf] * ":"
      tsv[key] = values

      tsv.add_field "Sample" do 
        sample
      end

      tsv.add_field "MAF" do 
        maf
      end

      tsv.add_field "BED" do 
        bed_type
      end

      if all.nil?
        all = tsv
      else
        all.merge!(tsv)
      end
    end
    all
  end

  dep :explore_values
  extension :png
  task :ggplot_image => :binary do
    require 'rbbt/util/R'
    require 'rbbt/util/R/plot'
    R::PNG.plot(self.tmp_path, step(:explore_values).load, <<-EOF)
names(data) <- make.names(names(data))
data$BED = factor(data$BED, levels=c('50', '100', '150', '200','capture'))
p = ggplot(data) + geom_boxplot(aes(x=factor(MAF), y=F.measure)) + facet_wrap(~BED)
print(p)
    EOF
    nil
  end

  input :vcf_file, :file, "VCF file to evaluate", nil, :required => true
  input :donor, :select, "Donor code", nil, :select_options => SAMPLES, :required => true
  dep :sample_vcfeval, :bed => :placeholder  do |jobname,options|
    vcf_file, donor = options.values_at :vcf_file, :donor
    donor ||= jobname
    options = options.merge :wgs => vcf_file, :bed => nil
    {:inputs => options, :jobname => donor}
  end
  task :evaluate => :tsv do
    step(:sample_vcfeval).load
  end

end

require 'tasks/genomic_mutations'

#require 'rbbt/knowledge_base/MODULE'
#require 'rbbt/entity/MODULE'

