all:
	echo "library(Rcpp); compileAttributes(\"./trajcomp\")" | R --no-save -q
	R CMD INSTALL trajcomp
build:
	R CMD build trajcomp
clean: build
	
