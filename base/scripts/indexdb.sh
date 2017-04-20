#!/bin/bash
mkdir index
for f in *.fasta; do 
 indexdb_rna --ref $f,./index/${f/.fasta/} &
done 

wait
