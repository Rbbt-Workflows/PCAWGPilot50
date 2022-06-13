module PCAWGPilot50

  task :merged_capture_bed_file => :text do
    CMD.cmd_log("cat '#{Rbbt.data.capture_beds.find}'/*.bed | bedtools sort | bedtools merge > #{self.tmp_path}")
    nil
  end

  dep :bed_file
  input :bases_to_pad, :integer, "Number of bases to add to regions", 200
  extension :bed
  task :padded_bed_file => :text do  |bases|
    reference = HTS.helpers[:reference_file].call 'hg38'
    reference = GATK.prepare_FASTA reference
    dict = reference.sub('.fa.gz', '.dict')

    sizes = {}
    Open.read(dict).split("\n").each do |line|
      next unless line =~ /SN:/
      name, size = line.split("\t").values_at 1, 2
      name.sub!("SN:", '')
      size.sub!("LN:", '')
      sizes[name] = size
    end

    TmpFile.with_file do |genome|
      Open.write(genome) do |f|
        sizes.each do |name,size|
          f.puts [name, size] * "\t"
        end
      end

      CMD.cmd(:bedtools, "slop -i #{step(:bed_file).path} -g #{genome} -b #{bases} | bedtools sort | bedtools merge  > #{self.tmp_path}")
    end
    nil
  end

  dep :padded_bed_file
  dep_task :slice_donor_bundle, HTSBenchmark, :slice_bundle, :regions_to_slice => :padded_bed_file do |donor,options|
    donor = donor + '-WGS-ARGO_aln'
    {:inputs => options, :jobname => donor}
  end

  dep :merged_capture_bed_file
  dep_task :merged_sliced_ref, HTSBenchmark, :sliceref, :bed_file => :merged_capture_bed_file, :reference => 'hg38', :do_vcf => true

  dep :merged_sliced_ref
  extension 'fa.gz'
  task :merged_sliced_ref_bundle => :binary do 

    log :reference, "Preparing reference"

    ref_dir = file('reference')
    slicedref = step(:merged_sliced_ref)
    slicedref.files_dir.glob("*/*.fa*").each do |file|
      Open.link file, ref_dir[File.basename(file)] 
    end

    slicedref.files_dir.glob("*/known_sites/*").each do |file|
      Open.link file, ref_dir.known_sites[File.basename(file)] 
    end

    HTSBenchmark.helpers[:prepare_FASTA].call(ref_dir.glob("*.fa.gz").first, ref_dir)

    ref_dir.known_sites.glob("*.vcf.gz").each do |vcf|
      GATK.prepare_VCF(vcf, ref_dir.known_sites)
    end

    Open.cp step(:merged_capture_bed_file).path, file('inputs')["regions.bed"]
  end
end
