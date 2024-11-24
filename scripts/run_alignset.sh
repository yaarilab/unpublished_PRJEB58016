#!/bin/bash

# Check if a file is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <input_fastq_file>"
    exit 1
fi

# Input FASTQ file
input_fastq=$1
output_fastq="${input_fastq%.fastq}_align-pass.fastq"

# Create a directory for logs if it doesn't exist
mkdir -p tmp_logs
mkdir -p tmp_logs2
mkdir -p tmp
# Extract unique UMIs
umis=$(grep '^@' $input_fastq | sed -n 's/.*UMI=\([A-Z]*-[0-9]\).*/\1/p' | sort | uniq)

# Count the number of unique UMIs
umi_count=$(echo "$umis" | wc -l)

# Print the number of unique UMIs
echo "Number of unique UMIs: $umi_count"

# Function to process a UMI in parallel
process_umi() {
    umi=$1
    tmp_fastq="tmp/tmp_${umi}.fastq"
    
    # Subset reads with the UMI and save to temporary file
    awk -v umi="$umi" 'BEGIN {FS = "\n"; OFS = "\n"} 
    {
        if ($0 ~ /^@/ && $0 ~ umi) {
            header = $0;
            getline seq;
            getline plus;
            getline qual;
            print header, seq, plus, qual;
        }
    }' $input_fastq > $tmp_fastq

    # Define muscle_exec_argv
    muscle_exec_argv="--exec /usr/local/bin/muscle"
    
    # Run AlignSets.py command
    AlignSets.py muscle -s $tmp_fastq --bf UMI $muscle_exec_argv --div --log "tmp_logs/tmp_${umi}.log" --nproc 1 --fail > "tmp_logs2/tmp_${umi}.log"

    # Cleanup the temporary FASTQ file after the process
    rm $tmp_fastq
}

# Function to monitor progress
monitor_progress() {
    while :; do
        log_count=$(ls -1 tmp_logs/ | wc -l)
        progress=$(echo "scale=2; ($log_count/$umi_count) * 100" | bc)
        echo "Progress: $log_count/$umi_count files processed ($progress%)"
        sleep 10
        # Exit if the process is complete
        if [ "$log_count" -ge "$umi_count" ]; then
            echo "All UMIs processed!"
            break
        fi
    done
}

# Start monitoring progress in the background
monitor_progress &

# Process each UMI in parallel
for umi in $umis; do
    process_umi "$umi" &
done

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
