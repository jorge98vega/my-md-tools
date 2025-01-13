for i in {1..49} # window
do
b=$(echo "scale=3; -1.550+($i*0.050)" | bc) # centro del potencial
# de -1.5 a 2.5

mkdir -p wham${i}
cp run_wham.sh wham${i}
cp MD_wham.sh wham${i}

cat << EOF > wham${i}/wham${i}_rst.dat
!!! RESTRAINED COORDINATES !!!
! 1 - REACTION COORDINATE - d(NZ1-HZ) - d(HZ-O)
&rst iat = 3503, 3506, 3506, 5944,
      rstwt = 1, -1,
      r1 = -500, r2 = ${b}, r3 = ${b}, r4 = 500,
      rk2 = 400, rk3 = 400,
      &end
!!! FREE COORDINATES !!!
! 2 - DISTANCE - d(NZ1-HZ)
&rst iat = 3503, 3506,
      r1 = -500, r2 = ${b}, r3 = ${b}, r4 = 500,
      rk2 = 0, rk3 = 0,
      &end
! 3 - DISTANCE - d(HZ-O)
&rst iat = 3506, 5944,
      r1 = -500, r2 = ${b}, r3 = ${b}, r4 = 500,
      rk2 = 0, rk3 = 0,
      &end
! 4 - DISTANCE - d(NZ1-O)
&rst iat = 3503, 5944,
      r1 = -500, r2 = ${b}, r3 = ${b}, r4 = 500,
      rk2 = 0, rk3 = 0,
      &end
! 5 - ANGLE - a(NZ1-HZ-O)
&rst iat = 3503, 3506, 5944,
      r1 = -500, r2 = ${b}, r3 = ${b}, r4 = 500,
      rk2 = 0, rk3 = 0,
      &end
! 6 - REACTION COORDINATE - d(O-H) - d(H-NZ2)
&rst iat = 5944, 5945, 5945, 2300
      rstwt = 1, -1,
      r1 = -500, r2 = ${b}, r3 = ${b}, r4 = 500,
      rk2 = 0, rk3 = 0,
      &end
! 7 - DISTANCE - d(O-H)
&rst iat = 5944, 5945,
      r1 = -500, r2 = ${b}, r3 = ${b}, r4 = 500,
      rk2 = 0, rk3 = 0,
      &end
! 8 - DISTANCE - d(H-NZ2)
&rst iat = 5945, 2300,
      r1 = -500, r2 = ${b}, r3 = ${b}, r4 = 500,
      rk2 = 0, rk3 = 0,
      &end
! 9 - DISTANCE - d(O-NZ2)
&rst iat = 5944, 2300,
      r1 = -500, r2 = ${b}, r3 = ${b}, r4 = 500,
      rk2 = 0, rk3 = 0,
      &end
! 10 - ANGLE - a(O-H-NZ2)
&rst iat = 5944, 5945, 2300,
      r1 = -500, r2 = ${b}, r3 = ${b}, r4 = 500,
      rk2 = 0, rk3 = 0,
      &end
EOF
done
