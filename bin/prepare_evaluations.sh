#!/bin/bash
for donor in `cat data/donors`; do echo $donor; drbbt hpc orchestrate workflow.rb donor_evaluation_batch -jn $donor -cl; done
