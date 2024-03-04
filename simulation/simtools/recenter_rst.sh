MOL=4t10s_run01
DIR=prod_0
RUN=prod

j=1 # Numero del frame
echo "trajin ${DIR}/${MOL}_${RUN}_${j}.rst" >> input
echo "center :1-320" >> input
echo "image" >> input
echo "trajout ${DIR}/${MOL}_${RUN}_${j}_recenter.rst" >> input

cpptraj ${DIR}/${MOL}.top input 

rm input
