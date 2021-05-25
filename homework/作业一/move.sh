#! /usr/bin/bash

for path in `find ./ -name *.txt`
do
    cp $path ~/nature/AGAC/data
done
cd ~/nature/AGAC/data
cat *.txt > total
cat total |tr -cs "[:alnum:]" "\n" |tr "[:upper:]" "[:lower:]" > all
cat all |sort |uniq > sort-uniq
wc -w all > result
wc -w sort-uniq >> result

