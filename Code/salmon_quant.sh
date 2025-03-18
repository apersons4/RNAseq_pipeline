#!/bin/bash
 data_dir="/media/sf_rnaseq/00_fastq"

 # Loop through the R1 files to get sample names
 for r1_file in ${data_dir}/*_R1_001.fastq.gz; do
 	base_name=$(basename ${r1_file} _R1_001.fastq.gz)
 	
 	r2_file="${data_dir}/${base_name}_R2_001.fastq.gz"
 	
 	echo "Processing sample ${r2_file}"
 	# Run the quant on both files
 	salmon quant -i elegans_index -l A \
 	-1 ${r1_file} \
        -2 ${r2_file}\
		 -p 8 --validateMappings --gcBias -o quants/${base_name}_quant 	
done
