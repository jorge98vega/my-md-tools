for i in {1..28} # step
do
b=$(echo "scale=3; 2.650-($i*0.050)" | bc)  # centro del potencial
# de 2.6 a 1.25

cat << EOF > rst_st${i}.dat
&rst iat = 92, 107,
      iresid = 0,
      ifvari = 0,
      r1 = 0, r2 = ${b}, r3 = ${b}, r4 = 500,
      rk2 = 500, rk3 = 500,
      &end
EOF
done
