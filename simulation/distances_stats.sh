MOL=4tubes_run01
DIR=prod_dist
RUN=prod

> ${DIR}/${MOL}_distances.dat
> ${DIR}/${MOL}_stats.dat

INI=1 # Primer archivo
FIN=5 # Ultimo archivo
COLS=( 7 8 1 9 10 2 11 12 3 13 14 4 15 16 5 17 18 6 ) # Columnas

i=$INI
while [[ $i -le $FIN ]]
do
    awk '{for(i=2; i<=NF; ++i) printf $i" "; print ""}' ${DIR}/${MOL}_${RUN}_${i}_rst.out >> ${DIR}/${MOL}_distances.dat

    i=$(( i + 1 ))
done

for j in "${COLS[@]}"
do
    awk -v N=$j '{ sum += $N; sumsq += ($N)^2} END { if (NR > 0) printf "%f %f \n", sum/NR, sqrt((sumsq - sum^2/NR)/NR) }' ${DIR}/${MOL}_distances.dat >> ${DIR}/${MOL}_stats.dat
done
