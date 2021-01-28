#!/bin/sh

gene=$1
surface=$2
line=`cat specific_site | awk -v var="$gene" -F"\t" '$1==var{print}'`
active_site=`echo $line | awk '{print $2}'`
aa_site=`echo $line | awk '{print $3}'`
python map.py ${gene} ref/pdb_entry variant/${gene}.aa $active_site $aa_site $surface

