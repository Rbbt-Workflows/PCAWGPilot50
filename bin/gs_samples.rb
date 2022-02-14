#!/usr/bin/env ruby

require 'rbbt-util'
require 'rbbt/util/simpleopt'
require 'set'

$0 = "rbbt #{$previous_commands*" "} #{ File.basename(__FILE__) }" if $previous_commands

options = SOPT.setup <<EOF

Find all samples for which we have data

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

gs = Rbbt.gold_standard.find(:lib)

files = TSV.setup({}, "Donor~caller#:type=:flat")

callers = Set.new
types = Set.new
gs.glob("*").each do |vc_dir|
  vc = File.basename(vc_dir)
  callers << vc
  vc_dir.glob("*.vcf").each do |file|
    name = File.basename(file, '.vcf')
    donor, _sep, type = name.partition("-")
    type = case type
           when "DS"
             "ARGO-aln"
           when "DS-orig"
             "rbbt"
           when "DS-rbbt"
             "rbbt.v1"
           end

    types << type

    files[donor] ||= []
    files[donor] << [type, vc] * ":"
  end
end

fields = types.to_a.collect{|t| callers.to_a.collect{|c| [t, c] * ":" }}.flatten
tsv = TSV.setup({}, :key_field => "Donor", :fields => fields, :type => :list)

files.each do |donor,list|
  list.each do |elem|
    type, vc = elem.split(":")
    f = [type,vc] * ":"
    tsv[donor] ||= [nil] * fields.length
    tsv[donor][fields.index(f)] = f
  end
end

puts tsv.to_s
