module PCAWGPilot50

  dep :slice_donor_bundle, :compute => :produce
  dep HTSBenchmark, :run_bundle, :bundle => :slice_donor_bundle, :input_type => "BAM"
  dep_task :evaluate_slice, self, :evaluate, :vcf_file => :run_bundle

end
