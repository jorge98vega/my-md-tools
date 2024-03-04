MOL=4t10s_run01
DIR=prod_0
RUN=prod

gunzip ${DIR}/${MOL}_${RUN}_*.coord.gz

j=1 #Contador de rmsd
while [[ $j -le 500 ]]
do
echo "trajin ${DIR}/${MOL}_${RUN}_${j}.coord " >> input
j=$(( j + 1 ))
done
echo "center :1-320" >> input
echo "image" >> input
echo "trajout ${DIR}/${MOL}_MD.nc netcdf" >> input

cpptraj ${DIR}/${MOL}.top input 

rm input
