INFILE="model_building/4t10s_runNoCl_prod_0_premove.rst"
INDECES=(6812) # The atom index if rst (starts at 0),
# the atom "serial" if pdb (starts at 1)
XYZ="newcoords.xyz" # The file with the new atom coordinates
OUTFILE="model_building/4t10s_runNoCl_prod_0.rst"
FORMAT="rst"

for INDEX in ${INDECES[@]}
do
if [ $FORMAT = "rst" ]
then
LINE=$((INDEX/2+3))
COLUMN=$((INDEX%2*3+1))

X=$(awk '{print $1}' $XYZ)
Y=$(awk '{print $2}' $XYZ)
Z=$(awk '{print $3}' $XYZ)

awk -v line=$LINE -v column=$COLUMN -v x=$X -v y=$Y -v z=$Z 'BEGIN {OFS="  "} NR==line {$column=x; $(column+1)=y; $(column+2)=z; $1="  "$1} 1' $INFILE > $OUTFILE
fi

if [ $FORMAT = "pdb" ]
then
LINE=$((INDEX+1))
COLUMN=7 # columna de la x, comprobar en el pdb, puede cambiar

X=$(awk '{print $1}' $XYZ)
Y=$(awk '{print $2}' $XYZ)
Z=$(awk '{print $3}' $XYZ)

awk -v line=$LINE -v column=$COLUMN -v x=$X -v y=$Y -v z=$Z 'NR==line {$column=x; $(column+1)=y; $(column+2)=z} 1' $INFILE > $OUTFILE
fi
done
