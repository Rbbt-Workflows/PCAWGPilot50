module PCAWGPilot50
  input :donor, :string, "Sample name", "DO50415"
  input :final, :boolean, "Use final list of mutations or the unfiltered", true
  dep_task :ds_genomic_mutations, Sequence, :genomic_mutations, :vcf => :placeholder do |jobname,options|
    donor = options[:donor]
    if options[:final]
      vcf = Rbbt.data.vcf["#{donor}-DS.vcf"].find
    else
      vcf = Rbbt.data.vcf_filtered["#{donor}-DS.vcf"].find
    end
    {:inputs => options.merge(:vcf_file => vcf), :task => :genomic_mutations}
  end

  input :donor, :string, "Sample name", "DO50415"
  input :final, :boolean, "Use final list of mutations or the unfiltered", true
  dep_task :wgs_genomic_mutations, Sequence, :genomic_mutations, :vcf => :placeholder do |jobname,options|
    donor = options[:donor]
    if options[:final]
      vcf = Rbbt.data.vcf["#{donor}-WGS.vcf"].find
    else
      vcf = Rbbt.data.vcf_filtered["#{donor}-WGS.vcf"].find
    end
    {:inputs => options.merge(:vcf_file => vcf), :task => :genomic_mutations}
  end

  dep :ds_genomic_mutations
  dep :wgs_genomic_mutations
  task :intersect_genomic_mutations => :array do
    ds = step(:ds_genomic_mutations).load
    wgs = step(:wgs_genomic_mutations).load
    ds & wgs
  end

  input :donor, :string, "Sample name", "DO50415"
  dep_task :lifted, Sequence, :lift_over, :positions => :placeholder, :source => "hg19", :target => "hg38" do |jobname,options|
    options[:positions] = PCAWG.genotypes[options[:donor]].find.list
    {:inputs => options}
  end

  dep :lifted
  dep :wgs_genomic_mutations
  task :intersect_lifted_wgs => :array do
    lifted = step(:lifted).load
    wgs = step(:wgs_genomic_mutations).load
    lifted & wgs
  end

  dep :lifted
  dep :ds_genomic_mutations
  task :intersect_lifted_ds => :array do
    lifted = step(:lifted).load
    ds = step(:ds_genomic_mutations).load
    lifted & ds
  end

  task :sizes => :tsv do
    tsv = TSV.setup({}, "Sample~Good_Mutations,All_Mutations#:type=:list")
    Rbbt.data.vcf.glob("*.vcf").sort.each do |file|
      donor = File.basename(file, '.vcf')
      tsv[donor] ||= []
      tsv[donor][0] = CMD.cmd("grep -v '#' #{file}|wc -l").read.strip.to_i
    end

    Rbbt.data.vcf_filtered.glob("*.vcf").sort.each do |file|
      donor = File.basename(file, '.vcf')
      tsv[donor] ||= []
      tsv[donor][1] = CMD.cmd("grep -v '#' #{file}|wc -l").read.strip.to_i
    end
    
    tsv
  end
end
