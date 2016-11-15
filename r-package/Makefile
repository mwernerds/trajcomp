all:
	echo "library(Rcpp); compileAttributes(\"./trajcomp\")" | R --no-save -q
	R CMD INSTALL trajcomp
touch:
	touch trajcomp/src/public_api.cpp
build:
	R CMD build trajcomp
clean: build
	
