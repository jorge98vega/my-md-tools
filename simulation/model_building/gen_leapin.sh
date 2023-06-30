INFILE=9tubes_10stack_oriented
OUTFILE=9t10s_run01
TUBES=9
RINGS=10
# ¡ATENCION! La secuencia es de 2 anillos, no de 1
SEQ="LYS PHD LYN PHD LYS PHD LYN PHD
LYN PHD LYS PHD LYN PHD LYS PHD"
# Iones por cada 2 anillos, se puede dejar vacío
IONS="TFA TFA TFA TFA"

cat << EOF > leap.in
source leaprc.protein.ff19SB
source leaprc.water.tip3p
verbosity 2
loadamberparams frcmod.ionsff99_tip3p
loadamberparams frcmod.tip3p
loadamberprep PHD.prep
loadamberprep TYD.prep
loadamberprep TFA.prepc
loadamberparams TFA.frcmod
x = loadpdbusingseq ${INFILE}.pdb {
EOF

for (( i=1; i<=$TUBES; i++ ))
do
for (( j=1; j<=$(( $RINGS / 2 )); j++ ))
do

echo $SEQ >> leap.in

done
done

for (( i=0; i<$TUBES; i++ ))
do
for (( j=1; j<=$(( $RINGS / 2 )); j++ ))
do 

echo $IONS >> leap.in

done
done

echo "}" >> leap.in

for (( i=0; i<$TUBES; i++ ))
do
for (( j=0; j<$RINGS; j++ ))
do

cat << EOF >> leap.in
bond x.$(( 1 + 8 * $j + 8 * $RINGS * $i )).N x.$(( 8 + 8*$j + 8 * $RINGS * $i )).C
EOF

done
done

cat << EOF >> leap.in
solvatebox x TIP3PBOX 15.0 2.0
addions x Cl- 180
addions x Na+ 0
charge x
savepdb x ${OUTFILE}.pdb
saveamberparm x ${OUTFILE}.top ${OUTFILE}.rst
EOF

# cat << EOF >> leap.in
# set x box { 40.0  40.0  40.0 }
# set default nocenter on
# charge x
# savepdb x ${OUTFILE}.pdb
# saveamberparm x ${OUTFILE}.top ${OUTFILE}.rst
# EOF
