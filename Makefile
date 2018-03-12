##### install the suite

.PHONY: all dependencies annotation gtf
all: init dependencies annotation
	@echo "Install complete"

init:
	@mkdir -p dependencies
	@mkdir -p annotation/gtf

## install dependencies
dependencies: dependencies/picard.jar
	
dependencies/picard.jar:
	wget https://github.com/broadinstitute/picard/releases/download/2.9.3/picard.jar -O $@

## install annotation
annotation: gtf

gtf: annotation/gtf/gencode.v19.annotation.gtf annotation/gtf/gencode.vM14.annotation.gtf

annotation/gtf/gencode.v19.annotation.gtf:
	wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz -O $@.gz
	gunzip $@.gz

annotation/gtf/gencode.vM14.annotation.gtf:
	wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M14/gencode.vM14.annotation.gtf.gz -O $@.gz
	gunzip $@.gz
