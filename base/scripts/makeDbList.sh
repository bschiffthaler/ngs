#!/bin/bash

dbFasta=$(find /usr/share/rRNA_databases -name "*.fasta")
dbLen=0
for i in $dbFasta; do let dbLen+=1; done

currDb=0
echo -n "SORTMERNA_DB="
for db in $dbFasta; do
    let currDb+=1;
    nam=$(basename $db);
    echo -n "$db,/usr/share/rRNA_databases/index/${nam/.fasta/}";
    [[ currDb -lt dbLen ]] && echo -n ":";
done
echo
