#!/usr/bin/env ruby

require 'rbbt-util'
require 'rbbt/util/simpleopt'
require 'rbbt/workflow'

$0 = "rbbt #{$previous_commands*" "} #{ File.basename(__FILE__) }" if $previous_commands

options = SOPT.setup <<EOF

Characterize FN

$ #{$0} [options] <filename.tsv|->

Use - to read from STDIN

-h--help Print this help

EOF
if options[:help]
  if defined? rbbt_usage
    rbbt_usage 
  else
    puts SOPT.doc
  end
  exit 0
end

wf = Workflow.require_workflow "PCAWGPilot48"

donor = ARGV.first || "DO50346"

normal = Rbbt.results["ARGO_mutect2/#{donor}-Sliced-WGS.vcf"].find
filters = Rbbt.results["ARGO_mutect2_filters/#{donor}-Sliced-WGS.vcf"].find

job = wf.job(:evaluate, nil, :donor => donor, :vcf_file => normal)
job.produce

positions = []
TSV.traverse job.step(:vcfeval).file('output/fn.vcf.gz'), :type => :list, :into => positions do |chr,values|
  position, *rest = values
  [chr, position]
end

causes = []
TSV.traverse filters, :type => :list, :into => causes do |chr,values|
  position, *rest = values
  next unless positions.include?([chr, position])
  rest[4]
end

Misc.counts(causes).sort_by{|c,count| count}.reverse.each do |cause,count|
  puts [cause, count] * "= "
end


