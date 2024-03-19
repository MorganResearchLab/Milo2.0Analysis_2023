#!/bin/bash
module load r/4.2.2

BASE_DIR="/nfs/research/marioni/mdmorgan/Milo2.0_Synthetic"
SRC=$(echo $BASE_DIR"/src")
DATA=$BASE_DIR"/data"
#pop_file=/nfs/team205/ed6/data/milo_benchmark/pop_sample_2_clean.txt
outdir=$BASE_DIR"/results"
logs=$BASE_DIR"/logs"

## Activate conda env
#conda activate milo_bm

# ## Run embryo
# cat $pop_file | \
# while read pop
#     do
#     for pop_enr in $(seq 0.75 0.1 0.95)
#         do
#         for seed in $(seq 43 1 45)
#             do
#             echo "Rscript ./make_bm_data.R --pop_enrichment $pop_enr /nfs/team205/ed6/data/milo_benchmark/embryo_data_bm.RDS $seed $pop" | bsub -o ${outdir}/milo_make_bm_data_${seed}.out -e ${outdir}/milo_make_bm_data_${seed}.err -G team283 -R"select[mem>3500] rusage[mem=3500]" -M3500
#             done
#         done
#     done

# ## Run synthetic
# data_id=linear3x
# for p in $(seq 1 1 7)
#     do 
#     for seed in 43 44 45
#         do 
#         for enr in $(seq 0.75 0.1 0.95)
#         do 
#             echo "Rscript ./make_bm_data.R --pop_enrichment $enr --make_batch_effect no --data_id $data_id /nfs/team205/ed6/data/milo_benchmark/${data_id}_data_bm.RDS ${seed} M${p}" | bsub -o ${outdir}/milo_make_bm_data_${data_id}_${seed}.out -e ${outdir}/milo_make_bm_data_${data_id}_${seed}.err -G team283 -R"select[mem>3500] rusage[mem=3500]" -M3500
#         done
#     done
# done

## Run synthetic
data_id=cluster3x
for p in $(seq 1)
    do 
    for seed in 43 44 45
        do 
        for enr in $(seq 0.7 0.1 0.9)
        do
	    JNAME="data_id_P"$p"_seed"$seed"_Enr"$enr

	    JOB="bsub -q standard  -M 3500 -R "rusage[mem=3500]" -J $JNAME -n1 -T 120  -e $logs/$NAME.err  -o $logs/$JNAME.out Rscript $SRC/make_bm_data_clusters.R --pop_enrichment $enr --population $p --batch_seed $seed  --make_batch_effect no --data_id $data_id --data_RDS $DATA/Clusters4_Ncells2000_20240122_discreteSim.RDS  --output $DATA/"

	    echo $JOB
	    eval $JOB
        done
    done
done

# ## Run synthetic
# data_id=branching3x
# for p in $(seq 1 1 10)
#     do 
#     for seed in 43 44 45
#         do 
#         for enr in $(seq 0.75 0.1 0.95)
#         do 
#             echo "Rscript ./make_bm_data.R --pop_enrichment $enr --make_batch_effect no --data_id $data_id /nfs/team205/ed6/data/milo_benchmark/${data_id}_data_bm.RDS ${seed} M${p}" | bsub -o ${outdir}/milo_make_bm_data_${data_id}_${seed}.out -e ${outdir}/milo_make_bm_data_${data_id}_${seed}.err -G team283 -R"select[mem>3500] rusage[mem=3500]" -M3500
#         done
#     done
# done


# if [ -f  ]; then
