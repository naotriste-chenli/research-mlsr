perc=$(awk -F ',' 'NR > 1 {print $7}' all_predictions.csv | datamash perc:95 1)
awk -F ',' 'NR > 1 {print $7 " " $1}' all_predictions.csv | awk -v x="$perc" '$1 > x {print $2}'