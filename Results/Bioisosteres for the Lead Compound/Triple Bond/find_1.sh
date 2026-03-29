#!/bin/env bash

perc=$(awk -F ',' '$7 != "QED" {print $7}' $1 | datamash perc:99 1)

df=$(awk -F ',' '$7 != "QED" && $1 != "smiles" {print $1 " " $7}' $1)

echo "$df" | awk -F ' ' -v x="$perc" '$2 >= x {print $1}'