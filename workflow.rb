require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/MODULE'

Workflow.require_workflow "HTS"
Workflow.require_workflow "Sequence"
Workflow.require_workflow "PCAWG"
module PCAWGPilot
  extend Workflow


  input :vcf, :file, "VCF file"
  input :bed, :file, "BED file"
  task :subset_by_bed do |vcf,bed|
    output = file('output/sample')
    Open.mkdir file('output')
    CMD.cmd("vcftools --vcf #{vcf} --bed #{bed}  --out #{output} --recode --keep-INFO-all")
    Open.mv file('output').glob("*.vcf").first, self.tmp_path
    nil
  end

  dep :subset_by_bed
  input :maf, :float, "Minimun maf", 0.1
  task :subset_by_af => :text do |maf|
    vcf = step(:subset_by_bed)

    TSV.traverse vcf, :type => :array, :into => :stream do |line|
      next line if line =~ /^#/
      parts = line.split("\t")

      format = parts[8].split(":")
      sample = parts[10].split(":")
      hash = Misc.zip2hash(format, sample)
      next if hash["AF"].to_f < maf
      line
    end
  end

  input :ds, :file
  input :wgs, :file
  dep :subset_by_af, :vcf => :ds, :jobname => "DS"
  dep :subset_by_bed, :vcf => :wgs, :jobname => "WGS"
  dep_task :vcfeval, HTS, :vcfeval, :input_vcf => :placeholder, :truth_vcf => :placeholder, :reference => 'hg38' do |jobname,options, dependencies|
    ds, wgs = dependencies.flatten.values_at 0, 1
    options[:input_vcf] = wgs
    options[:truth_vcf] = ds
    {:inputs => options, :jobname => jobname}
  end


end

require 'tasks/genomic_mutations'

#require 'rbbt/knowledge_base/MODULE'
#require 'rbbt/entity/MODULE'

