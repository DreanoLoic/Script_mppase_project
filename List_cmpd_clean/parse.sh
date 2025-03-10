#!/bin/bash
for i in {2..152}
do
   test=$(awk -F ',' ' FNR=='$i' {print $2}' ketamine_cmpd.csv);
   filename=$(grep -R -l "$test$" ./sdfs_divided_by_vendor/)
   filename="${filename//.\/sdfs_divided_by_vendor\//}"
   filename="${filename//.sdf/}"
   filename="${filename//.csv/}"
   filename="${filename//$'\n'/, }"
	echo "$test	$filename "  >> l_supplier.csv
done
