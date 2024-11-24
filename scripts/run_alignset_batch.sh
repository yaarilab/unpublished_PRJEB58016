#!/bin/bash

# Check if a file is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 <input_fastq_file> [batch_size]"
    exit 1
fi

# Input FASTQ file
input_fastq=$1
output_fastq="${input_fastq%.fastq}_align-pass.fastq"

# Create directories for logs and temporary files
mkdir -p tmp_logs
mkdir -p tmp_logs2
mkdir -p tmp

# Extract unique UMIs and group them into batches of XX
umis=$(grep '^@' $input_fastq | sed -n 's/.*UMI=\([A-Z]*-[0-9]\).*/\1/p' | sort | uniq)
umi_count=$(echo "$umis" | wc -l)

# Print the number of unique UMIs
echo "Number of unique UMIs: $umi_count"

# Function to process a batch of XX UMIs in parallel
process_batch() {
    batch_umis=$1
    batch_id=$2
    tmp_fastq="tmp/tmp_batch_${batch_id}.fastq"
    
    # Subset reads with the UMIs in the batch and save to a temporary file
    for umi in $batch_umis; do
        awk -v umi="$umi" 'BEGIN {FS = "\n"; OFS = "\n"} 
        {
            if ($0 ~ /^@/ && $0 ~ umi) {
                header = $0;
                getline seq;
                getline plus;
                getline qual;
                print header, seq, plus, qual;
            }
        }' $input_fastq >> $tmp_fastq
    done

    # Define muscle_exec_argv
    muscle_exec_argv="--exec /usr/local/bin/muscle"
    
    # Run AlignSets.py command
    AlignSets.py muscle -s $tmp_fastq --bf UMI $muscle_exec_argv --div --log "tmp_logs/tmp_batch_${batch_id}.log" --nproc 1 --fail > "tmp_logs2/tmp_batch_${batch_id}.log"

    # Cleanup the temporary FASTQ file after the process
    rm $tmp_fastq
}

# Function to monitor progress
monitor_progress() {
    while :; do
        log_count=$(ls -1 tmp_logs/ | wc -l)
        progress=$(echo "scale=2; ($log_count*100/($umi_count/$batch_size))" | bc)
        echo "Progress: $log_count/$(($umi_count/$batch_size)) batches processed ($progress%)"
        sleep 10
        # Exit if the process is complete
        if [ "$log_count" -ge "$(echo $umi_count/$batch_size | bc)" ]; then
            echo "All UMI batches processed!"
            break
        fi
    done
}

# Start monitoring progress in the background
monitor_progress &

# Split UMIs into batches of XX and process each batch in parallel
batch_id=0
batch_size=${2:-100}
batch_umis=""

for umi in $umis; do
    batch_umis+="$umi "
    if [ $(echo "$batch_umis" | wc -w) -eq $batch_size ]; then
        process_batch "$batch_umis" "$batch_id" &
        batch_umis=""
        batch_id=$((batch_id + 1))
    fi
done

# Process any remaining UMIs in the last batch
if [ ! -z "$batch_umis" ]; then
    process_batch "$batch_umis" "$batch_id" &
fi

# Wait for all background jobs to finish
wait

# Wait for the monitoring to finish
wait

echo "Concatenating all files into $output_fastq"
cat tmp/*.fastq > "$output_fastq"

# Cleanup temporary files
echo "Cleaning up temporary files"
rm -rf tmp tmp_logs tmp_logs2

echo "Alignment process complete. Output saved to $output_fastq"
