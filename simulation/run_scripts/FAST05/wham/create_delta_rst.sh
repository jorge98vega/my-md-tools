for i in {1..81} # step
do
b=$(echo "scale=3; -1.550+($i*0.050)" | bc) # centro del potencial
# de -1.5 a 2.5

cat << EOF > rst_st${i}.dat
&rst iat = 120, 107, 107, 92,
      ifvari = 0,
      rstwt = 1, -1,
      r1 = -500, r2 = ${b}, r3 = ${b}, r4 = 500,
      rk2 = 1000, rk3 = 1000,
      &end
EOF
done
