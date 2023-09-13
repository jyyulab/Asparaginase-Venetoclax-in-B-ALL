#!usr/bin/bash 

# Coded by Xin Huang(xhuang@stjude.org)

# python: 3.6.1
# SJARACNe: https://github.com/jyyulab/SJARACNe

indir = ./SJARACNe
cd $indir
for i in $(ls -d */ | cut -f1 -d'/'); do
echo sjaracne lsf -j $indir/config_cwlexec.json -e $indir/${i}/*.exp -g $indir/${i}/tf/*.txt -n 100 -o $indir/${i}/tf -pc 0.01;
done