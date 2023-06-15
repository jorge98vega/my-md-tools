rm -r umbrella_files
mkdir umbrella_files

for i in {1..28..1} # step
do
  b=$(echo "scale=3; 2.650-($i*0.050)" | bc) # centro del potencial
  awk '{print $1,$2}' rst_md${i}.out | tail -1500 > umbrella_files/mc_${i}.out
  echo mc_${i}.out ${b} 1000 >> umbrella_files/meta.dat # el doble que en rst
done

cd umbrella_files
/home/jorge/bin/wham/wham/wham 1.25 2.6 50 0.0001 300 0 meta.dat wham.out
# hist_min hist_max num_bins tolerance temperature periodic metafile outfile
cd ..
