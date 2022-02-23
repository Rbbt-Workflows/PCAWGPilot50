module PCAWGPilot50

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
end
