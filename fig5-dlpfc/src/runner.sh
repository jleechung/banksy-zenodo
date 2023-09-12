
# 11 seeds
all_seed=(0 1 2 10 41 42 100 123 1000 1234 2020)

for curr_seed in "${all_seed[@]}"; do
    Rscript banksy.R --sampleid 1 --seed $curr_seed &
    Rscript banksy.R --sampleid 2 --seed $curr_seed &
    Rscript banksy.R --sampleid 3 --seed $curr_seed &
    Rscript banksy.R --sampleid 4 --seed $curr_seed &
    Rscript banksy.R --sampleid 5 --seed $curr_seed &
    Rscript banksy.R --sampleid 6 --seed $curr_seed &
    Rscript banksy.R --sampleid 7 --seed $curr_seed &
    Rscript banksy.R --sampleid 8 --seed $curr_seed &
    Rscript banksy.R --sampleid 9 --seed $curr_seed &
    Rscript banksy.R --sampleid 10 --seed $curr_seed &
    Rscript banksy.R --sampleid 11 --seed $curr_seed &
    Rscript banksy.R --sampleid 12 --seed $curr_seed &
    wait
done

