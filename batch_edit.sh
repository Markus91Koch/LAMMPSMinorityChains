#!/bin/sh
# this script takes a brush in lammps data format (atom style bond)
# and extends one chain at the center (via python script)
# plus edits the necessary header entries

if [ -z "$1" ]
  then
    echo "Please provide the filename of the original lammps file as the 1st input parameter."
    exit 1
fi

# read in and split filename
extension="${1##*.}"
filename="${1%.*}"
start_fn="${1%_*}"
end_fn="${1##*_}"

# get certain line numbers of original lammps file for splitting
# line number where "Atoms" and "Bonds" keyword is found
na=$(grep -n Atoms $1 |cut -f1 -d:)
nb=$(grep -n Bonds $1 |cut -f1 -d:)
nl=$(wc -l $1 | awk '{print $1}')       # total length of file

# split orig. file into header, footer and atoms section
head -n+$((na+1)) $1 > header.txt
tail -n+$((nb+2)) $1 > footer.txt
tail -n+$((nb-1)) $1 | head -3 > mid.txt
tail -n+$((na+2)) $1 | head -$((nb-na-3)) > atoms.txt


for k in 192 256 320
do
    # apply python script on brush to extend central chain
    # with special end groups (separate atom type)
    echo "$k"
    python minority.py "$k"

    # old and new number of atoms
    oldN=$(grep atoms header.txt | cut -d ' ' -f 1)
    newN=$(head -1 newinfo.txt)

    # old and new number of bonds
    oldB=$(grep bonds header.txt | cut -d ' ' -f 1)
    newB=$(tail -1 newinfo.txt)

    # edit header of new data file
    sed -e "s/"$oldN"/"$newN"/g" header.txt > buffer.txt
    sed -e "s/"$oldB"/"$newB"/g" buffer.txt > newheader.txt

    # generate new and unique filename
    newname="$start_fn""_"$k"_""$end_fn"

    # stitch together all file sections to creat new data file
    cat newheader.txt newatoms.txt mid.txt newfooter.txt > "$newname"


done

# cleanup
rm newheader.txt newatoms.txt newfooter.txt buffer.txt newinfo.txt
rm header.txt footer.txt atoms.txt mid.txt

