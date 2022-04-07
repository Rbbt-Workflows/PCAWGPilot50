#!/bin/bash

for donor in `cat data/donors`; do echo $donor; drbbt task workflow.rb ds_consensus_vcf -jn $donor --mutation_type both --min_callers=3 --ds_system realn -O gold_standard/consensus_min3/$donor-DS-realn.vcf & done
for donor in `cat data/donors`; do echo $donor; drbbt task workflow.rb ds_consensus_vcf -jn $donor --mutation_type both --min_callers=2 --ds_system realn -O gold_standard/consensus_min2/$donor-DS-realn.vcf & done
for donor in `cat data/donors`; do echo $donor; drbbt task workflow.rb ds_consensus_vcf -jn $donor --mutation_type both --min_callers=1 --ds_system realn -O gold_standard/consensus_min1/$donor-DS-realn.vcf & done

for donor in `cat data/donors`; do echo $donor; drbbt task workflow.rb wgs_consensus_vcf -jn $donor --mutation_type both --min_callers=3 --wgs_system sliced -O result/consensus_min3/$donor-WGS-Sliced.vcf & done
for donor in `cat data/donors`; do echo $donor; drbbt task workflow.rb wgs_consensus_vcf -jn $donor --mutation_type both --min_callers=2 --wgs_system sliced -O result/consensus_min2/$donor-WGS-Sliced.vcf & done
for donor in `cat data/donors`; do echo $donor; drbbt task workflow.rb wgs_consensus_vcf -jn $donor --mutation_type both --min_callers=1 --wgs_system sliced -O result/consensus_min1/$donor-WGS-Sliced.vcf & done
