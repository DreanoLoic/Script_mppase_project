#!/bin/bash
for i in {1..5}
do
   awk -F ',' ' FNR==$i {print $1}' ketamine_cmpd.csv
done
