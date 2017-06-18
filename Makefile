##### install the suite

.PHONY: all dependencies
all: dependencies

## install dependencies
dependencies: dependencies/picard.jar
	
dependencies/picard.jar:
	wget https://github.com/broadinstitute/picard/releases/download/2.9.3/picard.jar -O $@
