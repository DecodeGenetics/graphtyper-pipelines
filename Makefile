test: clean
	samtools faidx test/reference.fa
	samtools index test/SAMP1.bam
	samtools index test/SAMP2.bam
	samtools index test/SAMP3.bam
	samtools index test/SAMP4.bam
	bash make_graphtyper_pipeline.sh test/bamlist test_config.sh | bash

clean:
	rm -rf test/results


.PHONY: test
