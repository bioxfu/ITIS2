source activate itis2
conda env export > doc/environment.yml

if [ ! -d fastq ]; then
	mkdir -p blast fastq fastq genome sam table TE
fi
