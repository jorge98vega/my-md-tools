MOL=4t10s_run01

echo "trajin prod_0/${MOL}_MD.nc" >> input
echo "trajin prod_1/${MOL}_MD.nc" >> input
echo "center :1-320" >> input
echo "image" >> input
echo "trajout ${MOL}_MD_merged.nc netcdf" >> input

cpptraj ${MOL}.top input 

rm input
