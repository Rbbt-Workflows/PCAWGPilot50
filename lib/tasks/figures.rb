module PCAWGPilot50
  dep :donor_vcfeval, :compute => :bootstrap, :maf => :placeholder, :bed_type => :placeholder do |jobname,options,dependencies|
    SAMPLES.collect do |donor|
      [0.1, 0.05, 0.01, 0.005, 0.001].collect do |maf|
        %w(capture 50 100 150 200).collect do |bed_type|
          {:inputs => options.merge(:maf => maf, :bed_type => bed_type), :jobname => donor}
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
      donor = dep.clean_name

      values = tsv[tsv.keys.first]

      tsv.keys.each{|k| tsv.delete k}
      key = [donor, bed_type, maf] * ":"
      tsv[key] = values

      tsv.add_field "Sample" do 
        donor
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
  input :statistic, :select, "Statistic to display", "F.measure", :select_options => %w(F.measure Precision Sensitivity)
  extension :png
  task :explore_values_figure => :binary do |statistic|
    require 'rbbt/util/R'
    require 'rbbt/util/R/plot'
    R::PNG.plot(self.tmp_path, step(:explore_values).load, <<-EOF)
names(data) <- make.names(names(data))
data$BED = factor(data$BED, levels=c('50', '100', '150', '200','capture'))
p = ggplot(data) + geom_boxplot(aes(x=factor(MAF), y=#{statistic})) + facet_wrap(~BED)
print(p)
    EOF
    nil
  end

  #{{{ CALLER EXPLORATION
  dep :donor_vcfeval, :compute => :bootstrap, :ds_caller => :placeholder, :wgs_caller => :placeholder, :donor => :placeholder do |jobname, options|
    SAMPLES.collect do |donor|
      %w(mutect2 strelka muse).collect do |ds_caller|
        %w(mutect2 strelka muse).collect do |wgs_caller|
          job_options = options.merge({:donor => donor, :ds_caller => ds_caller, :wgs_caller => wgs_caller})
          {:inputs => job_options, :jobname => donor}
        end
      end
    end.flatten
  end
  task :caller_exploration => :tsv do
    all = nil
    dependencies.each do |dep|
      tsv = dep.load
      donor = dep.clean_name
      ds_caller = dep.recursive_inputs[:ds_caller]
      wgs_caller = dep.recursive_inputs[:wgs_caller]

      values = tsv[tsv.keys.first]

      tsv.keys.each{|k| tsv.delete k}
      key = [donor, ds_caller, wgs_caller] * ":"
      tsv[key] = values

      tsv.add_field "Sample" do 
        donor
      end

      tsv.add_field "DS Caller" do 
        ds_caller
      end

      tsv.add_field "WGS Caller" do 
        wgs_caller
      end

      if all.nil?
        all = tsv
      else
        all.merge!(tsv)
      end
    end
    all
  end

  dep :caller_exploration
  input :statistic, :select, "Statistic to display", "F.measure", :select_options => %w(F.measure Precision Sensitivity)
  extension :png
  task :caller_exploration_figure => :binary do |statistic|
    require 'rbbt/util/R'
    require 'rbbt/util/R/plot'
    R::PNG.plot(self.tmp_path, step(:caller_exploration).load, <<-EOF)
names(data) <- make.names(names(data))
p = ggplot(data) + geom_boxplot(aes(x=factor(WGS.Caller), y=#{statistic})) + facet_wrap(~DS.Caller)
print(p)
    EOF
    nil
  end


  input :donor, :select, "Donor to compare to", nil, :select_options => SAMPLES, :jobname => true
  dep :ds_vcf, :ds_caller => "mutect2"
  dep :ds_vcf, :ds_caller => "strelka"
  dep :ds_vcf, :ds_caller => "muse"
  extension :png
  task :donor_ds_caller_venn => :binary do
    mutations = TSV.setup({}, "Mutation~mutect2,strelka,muse#:type=:list")
    dependencies.each do |dep|
      dep.join
      Sequence.job(:genomic_mutations, nil, :vcf_file => dep.path).run.each do |mut|
        mutations[mut] ||= NamedArray.setup([false, false, false], mutations.fields)
        mutations[mut][dep.recursive_inputs[:ds_caller]] = true
      end
    end
    R::PNG.plot(self.tmp_path, mutations, <<-EOF, 700, 700)
rbbt.require('gridExtra')
plot = rbbt.plot.venn(data)
plot = grid.arrange(gTree(children=plot), top="#{clean_name}")
print(plot)
    EOF
    nil
  end

  input :donor, :select, "Donor to compare to", nil, :select_options => SAMPLES, :jobname => true
  dep :ds_vcf, :ds_caller => "mutect2"
  dep :ds_vcf, :ds_caller => "strelka"
  dep :ds_vcf, :ds_caller => "muse"
  extension :png
  task :donor_ds_caller_upset => :binary do
    mutations = TSV.setup({}, "Mutation~mutect2,strelka,muse#:type=:list")
    dependencies.each do |dep|
      dep.join
      Sequence.job(:genomic_mutations, nil, :vcf_file => dep.path).run.each do |mut|
        mutations[mut] ||= NamedArray.setup([false, false, false], mutations.fields)
        mutations[mut][dep.recursive_inputs[:ds_caller]] = true
      end
    end
    R::PNG.plot(self.tmp_path, mutations, <<-EOF, 1100, 700)
plot = rbbt.plot.upset(data, text.scale=2)
print(plot)
    EOF
    nil
  end

  input :donor, :select, "Donor to compare to", nil, :select_options => SAMPLES, :jobname => true
  dep :wgs_vcf, :wgs_caller => "mutect2"
  dep :wgs_vcf, :wgs_caller => "strelka"
  dep :wgs_vcf, :wgs_caller => "muse"
  extension :png
  task :donor_wgs_caller_upset => :binary do
    mutations = TSV.setup({}, "Mutation~mutect2,strelka,muse#:type=:list")
    dependencies.each do |dep|
      dep.join
      Sequence.job(:genomic_mutations, nil, :vcf_file => dep.path).run.each do |mut|
        mutations[mut] ||= NamedArray.setup([false, false, false], mutations.fields)
        mutations[mut][dep.recursive_inputs[:wgs_caller]] = true
      end
    end
    R::PNG.plot(self.tmp_path, mutations, <<-EOF, 1100, 700)
plot = rbbt.plot.upset(data, text.scale=2)
print(plot)
    EOF
    nil
  end

  dep :donor_vcfeval, :compute => :canfail do |donor,options|
    jobs = []
    PCAWGPilot50::DS_CALLERS.each do |ds_caller|
      %w(realn ARGO-aln).each do |ds_system|
        PCAWGPilot50::WGS_CALLERS.each do |wgs_caller|
          %w(normal sliced).each do |wgs_system|
            %w(SNV indel both).each do |type|
    
              begin
                jobs << PCAWGPilot50.job(:donor_vcfeval, donor, 
                                         :wgs_caller => wgs_caller,
                                         :wgs_system => wgs_system,
                                         :ds_caller => ds_caller,
                                         :ds_system => ds_system,
                                         :mutation_type => type)
              rescue
              end
            end
          end
        end
      end
    end
    jobs
  end
  task :donor_evaluation_batch => :tsv do
    tsv = TSV.setup({}, "ID~Donor,Mutation type,DS caller, DS system,WGS caller,WGS system,TP,FP,FN,Precision,Sensitivity,F-score,vcfeval")
    donor = self.clean_name
    dependencies.each do |dep|
      next if dep.error? || dep.aborted?

      res = dep.load

      next if res.empty?

      values = [donor] + dep.recursive_inputs.values_at(:mutation_type, :ds_caller, :ds_system, :wgs_caller, :wgs_system)
      values += res[res.keys.first].values_at "True-pos-call", "False-pos", "False-neg", "Precision", "Sensitivity", "F-measure"
      values << dep.step(:vcfeval).path
      id = Misc.digest(values.inspect)
      tsv[id] = values

    end
    tsv
  end

  dep :donor_evaluation_batch do |jobname,options|
    SAMPLES.collect do |donor|
      {:inputs => options, :jobname => donor}
    end
  end
  task :evaluation_batch => :tsv do
    dependencies.inject(nil) do |acc,dep|
      tsv = dep.load
      if acc.nil?
        acc = tsv
      else
        acc.merge!(tsv)
        acc
      end
    end
  end

  input :evaluations, :tsv, "Evaluation values"
  extension :svg
  task :sensitivity_recall => :text do |evaluations|
    evaluations = TSV.open(evaluations) unless TSV === evaluations
    R::SVG.ggplot evaluations, <<-EOF
ggplot(data) + geom_point(aes(x=Sensitivity, y=Precision))
    EOF
  end

end
