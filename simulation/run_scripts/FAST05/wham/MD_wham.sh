#######################################################################
# DEFINE LOCALS
set -e

module load ips/2019

INPUT=$1
WINDOW=$2
ITER=$3
NSTEPS=$4
JOB=${INPUT}_wham${2}_${ITER}
PREV=$((ITER-1))

#######################################################################
# RUN SIMULATION

cat << EOF > amberwham.in
QM/MM
 &cntrl
 ifqnt = 1,
 nmropt = 1, ! (Leer fichero de restraints)
 ntx = 5, irest = 1, ntrx = 1, ! (Formato del input)
 ntxo = 2, iwrap = 0, ntpr = 1, ntwx = 1, ntwr = 1, ! (Formato y frecuencia del output)
 ntr = 0, ! (comentado) restraintmask = '@339', restraint_wt = 20, ! (Restraints sobre átomos)
 ntp = 2, ntf = 2, ntc = 2, ntb = 2, cut = 10, ! (Presión, SHAKE, PBCs)
 imin = 0, ! (Flags de minimización)
 nstlim = $NSTEPS, dt = 0.0005, ! (Flags de dinámica molecular) dt en ps -> 0.5 fs
 temp0 = 300, tempi = 300, ntt = 3, ig=-1, gamma_ln = 5, ! (Control de la temperatura)
 &end
 &wt type = 'DUMPFREQ', istep1 = 1, /
 &wt type = 'END' &end
 DISANG = wham${WINDOW}_rst.dat
 LISTOUT = wham${WINDOW}_rst_${ITER}.lis
 DUMPAVE = wham${WINDOW}_rst_${ITER}.dump
/
 &qmmm
 qmmask= '@761,762,763,764,765,766,767,3334,3335,3336,3337,3338,3339,3340,927,928,929,930,931,932,933,3500,3501,3502,3503,3504,3505,3506,2048,2049,2050,2051,2052,2053,4787,4788,4789,4790,4791,4792,2297,2298,2299,2300,2301,2302,4870,4871,4872,4873,4874,4875,5788,5789,5790,5842,5843,5844,6214,6215,6216,6226,6227,6228,6394,6395,6396,6439,6440,6441,6409,6410,6411,5944,5945,5946,6427,6428,6429,5376,5377,5378,5379,5380,5381,5382,5390,5391,5392,5393,5394,5395,5396,5593,5594,5595,5596,5597,5598,5599,5607,5608,5609,5610,5611,5612,5613'
 qm_theory='extern'
 qm_ewald= 0,
 qmshake= 0,
 qmcharge= 0,
 writepdb= 1,
/
 &fb
 basis= '/home/jorge/Fireball/Fdatas/Fdata_HCNOSF_47/',
 idipole = 1,
 idftd3  = 1
/
EOF

$AMBERHOME/bin/sander -O -i amberwham.in \
    -o ${JOB}.out \
    -p ${INPUT}.top \
    -c ${INPUT}_wham${WINDOW}_${PREV}.rst \
    -x ${JOB}.nc \
    -r ${JOB}.rst > ${JOB}.out

if [[ $ITER -eq 1 ]]; then
    if [[ $5 -eq 0 ]]; then
	NEXT_WINDOW=$((WINDOW-1))
        cp ${JOB}.rst ../wham${NEXT_WINDOW}/${INPUT}_wham${NEXT_WINDOW}_0.rst
	NEXT_WINDOW=$((WINDOW+1))
	cp ${JOB}.rst ../wham${NEXT_WINDOW}/${INPUT}_wham${NEXT_WINDOW}_0.rst
    elif [[ $5 -lt 0 ]]; then
	NEXT_WINDOW=$((WINDOW-1))
	cp ${JOB}.rst ../wham${NEXT_WINDOW}/${INPUT}_wham${NEXT_WINDOW}_0.rst
    elif [[ $5 -gt 0 ]]; then
	NEXT_WINDOW=$((WINDOW+1))
	cp ${JOB}.rst ../wham${NEXT_WINDOW}/${INPUT}_wham${NEXT_WINDOW}_0.rst
fi
