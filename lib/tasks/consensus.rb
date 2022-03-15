module PCAWGPilot50
  WGS_CONSENSUS_CALLERS =<<-EOF.split("\n")
ARGO_mutect2
SAGE
lofreq
somatic_sniper
varscan
  EOF

  DS_CONSENSUS_CALLERS =<<-EOF.split("\n")
muse
mutect2_pon_4.2.5
strelka
  EOF


  input :consensus_callers, :array, "Consensus callers to consider", PCAWGPilot50::WGS_CONSENSUS_CALLERS
  dep :wgs_vcf, :wgs_caller => :placeholder, :compute => :canfail do |jobname,options|
    options[:consensus_callers].collect do |ccaller|
      {:inputs => options.merge(:wgs_caller => ccaller), :jobname => jobname }
    end
  end
  extension :vcf
  task :wgs_combined_vcf => :text do
    files = {}
    dependencies.each do |dep|
      vcaller = dep.recursive_inputs[:wgs_caller]
      files[vcaller] = dep.path
    end
    HTS.combine_caller_vcfs(files)
  end

  dep :wgs_combined_vcf
  input :min_callers, :integer, "Min number of callers to pass variant", 2
  task :wgs_consensus_vcf => :text do |min_callers|
    TSV.traverse step(:wgs_combined_vcf), :type => :array, :into => :stream do |line|
      next line if line =~ /^#/
      parts = line.split("\t")
      filter = parts[6]

      callers = parts[6].split(";")
      callers = callers.select{|f| f.split("--").last == "PASS"}
      callers = callers.collect{|f| f.split("--").first }
      num = callers.uniq.length
      next unless num >= min_callers
      parts[6] = "PASS"
      parts * "\t"
    end
  end



  input :consensus_callers, :array, "Consensus callers to consider", PCAWGPilot50::DS_CONSENSUS_CALLERS
  dep :ds_vcf, :ds_caller => :placeholder, :compute => :canfail do |jobname,options|
    options[:consensus_callers].collect do |ccaller|
      begin
        PCAWGPilot50.job(:ds_vcf, jobname, options.merge(:ds_caller => ccaller))
      rescue
        next nil
      end
    end.compact
  end
  extension :vcf
  task :ds_combined_vcf => :text do
    files = {}
    dependencies.each do |dep|
      vcaller = dep.recursive_inputs[:ds_caller]
      files[vcaller] = dep.path
    end
    HTS.combine_caller_vcfs(files)
  end

  dep :ds_combined_vcf
  input :min_callers, :integer, "Min number of callers to pass variant", 2
  task :ds_consensus_vcf => :text do |min_callers|
    TSV.traverse step(:ds_combined_vcf), :type => :array, :into => :stream do |line|
      next line if line =~ /^#/
      parts = line.split("\t")
      filter = parts[6]

      callers = parts[6].split(";")
      callers = callers.select{|f| f.split("--").last == "PASS"}
      callers = callers.collect{|f| f.split("--").first }
      num = callers.uniq.length
      next unless num >= min_callers
      parts[6] = "PASS"
      parts * "\t"
    end
  end
end
