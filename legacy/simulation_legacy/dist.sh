DIR=estab_3
MOL=4tubes_run01
RUN=mdestab

rm ${DIR}/${MOL}_distance.dat

j=1  # Contador estructura final
while [[ $j -le 500 ]]
do
awk '{for(i=2; i<=NF; ++i) printf $i" "; print ""}' ${DIR}/${MOL}_${RUN}_${j}_rst.out >> ${DIR}/${MOL}_distance.dat
j=$(( j + 1 ))
done
