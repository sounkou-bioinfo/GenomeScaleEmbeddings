# h/t to @jimhester and @yihui for this parse block:
# https://github.com/yihui/knitr/blob/dc5ead7bcfc0ebd2789fe99c527c7d91afb3de4a/Makefile#L1-L4
# Note the portability change as suggested in the manual:
# https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Writing-portable-packages
PKGNAME = `sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION`
PKGVERS = `sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION`


.PHONY: all build document check install install_deps clean local_embeddings_duck local_houba
all: check

build: install_deps document
	R CMD build .

document:
	Rscript -e 'if (!requireNamespace("roxygen2")) install.packages("roxygen2")' \
	        -e 'roxygen2::roxygenise()'
check: build
	R CMD check --no-manual $(PKGNAME)_$(PKGVERS).tar.gz

install_deps:
	Rscript \
	-e 'if (!requireNamespace("remotes")) install.packages("remotes")' \
	-e 'remotes::install_deps(dependencies = TRUE)'

install: build
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

local_embeddings.duckdb: install
	R -e 'GenomeScaleEmbeddings::CopyParquetToDuckDB()'

local_embeddings.houba: local_embeddings.duckdb
	R -e 'GenomeScaleEmbeddings::writeEmbeddingsHoubaFromDuckDB()'

clean:
	@rm -rf $(PKGNAME)_$(PKGVERS).tar.gz $(PKGNAME).Rcheck
