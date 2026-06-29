PACKAGE := $(shell grep '^Package:' DESCRIPTION | sed -E 's/^Package:[[:space:]]+//')
RSCRIPT = Rscript --no-init-file

all: compile

compile: attributes
	${RSCRIPT} -e 'pkgbuild::compile_dll()'
	make roxygen

test:
	${RSCRIPT} -e 'library(methods); devtools::test()'

attributes:
	${RSCRIPT} -e "Rcpp::compileAttributes()"

roxygen:
	@mkdir -p man
	${RSCRIPT} -e "library(methods); devtools::document()"

install:
	R CMD INSTALL .

build:
	R CMD build --no-build-vignettes .

check: build
	R CMD check --no-build-vignettes --no-manual `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -f `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -rf ${PACKAGE}.Rcheck

clean:
	rm -f src/*.o src/*.so

.PHONY: all compile attributes roxygen test install build check clean
