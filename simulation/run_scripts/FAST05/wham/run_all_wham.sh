#/bin/bash

# Como usar el script - windows de 1 a 49, empezando por la 25:
# ./run_all_wham.sh 4t8sXwL_run40w 25 50 0

INPUT=$1
window=$2
limit=$3
step=$4
OLD_JOB_ID=$5
wd=$(pwd)

if [ $step -eq 0 ]; then
    cp ${wd}/${INPUT}.top ${wd}/wham${window}/
    cp ${wd}/${INPUT}_wham_0.rst ${wd}/wham${window}/${INPUT}_wham${window}_0.rst
    cd ${wd}/wham${window}/
    echo "Window $window - Step $step"
    JOB_ID=$(sbatch -J "prep${window}" run_wham.sh $INPUT $window 1 2000 $step | awk '{print $4}')
    sbatch -J "wham${window}" -d afterany:$((JOB_ID)) run_wham.sh $INPUT $window 2 5000
    echo "Job $JOB_ID"
    cd $wd

    new_window=$((window-1))
    new_step=$((step-1))
    ./run_all_wham.sh $new_window $new_step $((JOB_ID))

    new_window=$((window+1))
    new_step=$((step+1))
    ./run_all_wham.sh $new_window $new_step $((JOB_ID))
else; then
    cp ${wd}/${INPUT}.top ${wd}/wham${window}/
    cd ${wd}/wham${window}/
    echo "Window $window - Step $step"
    JOB_ID=$(sbatch -J "prep${window}" -d afterok:${OLD_JOB_ID} run_wham.sh $INPUT $window 1 2000 $step | awk '{print $4}')
    sbatch -J "wham${window}" -d afterany:$((JOB_ID)) run_wham.sh $INPUT $window 2 5000
    echo "Job $JOB_ID"
    cd $wd

    if [[ $step -lt 0 ]] && [[ $window -gt 0 ]]; then
    	new_window=$((window-1))
    	new_step=$((step-1))
        ./run_all_wham.sh $new_window $limit $new_step $((JOB_ID))
    elif [[ $step -gt 0 ]] && [[ $window -lt $limit ]]; then
	new_window=$((window+1))
        new_step=$((step+1))
        ./run_all_wham.sh $new_window $limit $new_step $((JOB_ID))
    fi
fi
