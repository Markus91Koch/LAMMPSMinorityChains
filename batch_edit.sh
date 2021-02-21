#!/bin/sh

if [ -z "$1" ]
  then
    echo "Please provide the filename of the original lammps file as the 3rd input parameter."
    exit 1
fi

extension="${1##*.}"
filename="${1%.*}"
start_fn="${1%_*}"
end_fn="${1##*_}"

# get certain line numbers of original lammps file
# line number where "Atoms" keyword is found
na=$(grep -n Atoms $1 |cut -f1 -d:)
nb=$(grep -n Bonds $1 |cut -f1 -d:)
nl=$(wc -l $1 | awk '{print $1}')

# split orig. file into header, footer and atoms section
head -n+$((na+1)) $1 > header.txt
tail -n+$((nb+2)) $1 > footer.txt
tail -n+$((nb-1)) $1 | head -3 > mid.txt
tail -n+$((na+2)) $1 | head -$((nb-na-3)) > atoms.txt

#sort -t, -nk1 atoms.txt -o atoms.txt
#sort -t, -nk1 velo.txt -o velo.txt

#cat header.txt atoms.txt mid.txt velo.txt footer.txt > final.data

#for k in 64 96 128 160 192 224 256 288 320 352 384 #64 80 96 112 128 144 160 176 192 
#for k in 96 128 192 256
for k in 192 256 320
do
    echo "$k"
    python minority.py "$k"

    #exit

    oldN=$(grep atoms header.txt | cut -d ' ' -f 1)
    newN=$(head -1 newinfo.txt)

    oldB=$(grep bonds header.txt | cut -d ' ' -f 1)
    newB=$(tail -1 newinfo.txt)

    sed -e "s/"$oldN"/"$newN"/g" header.txt > buffer.txt
    sed -e "s/"$oldB"/"$newB"/g" buffer.txt > newheader.txt

    newname="$start_fn""_"$k"_""$end_fn"

    cat newheader.txt newatoms.txt mid.txt newfooter.txt > "$newname" #"$filename"_sorted.data
    #rm newheader.txt newatoms.txt mid.txt newfooter.txt buffer.txt


done
rm newheader.txt newatoms.txt newfooter.txt buffer.txt newinfo.txt
rm header.txt footer.txt atoms.txt mid.txt

