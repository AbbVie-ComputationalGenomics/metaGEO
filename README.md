# CompGenomics-metaGEO

## R package requirements
It is recommended to install these packages via Bioconductor's biocLite() method:
* shiny
* shinyjs
* shinyBS
* lattice
* logging
* hgu133plus2.db
* org.Hs.eg.db
* knitr
* gdata
* dbConnect
* RSQLite
* GEOmetadb
* xtable
* rmarkdown
* devtools

The metageo object building process (see the [vignettes](vignettes/)) requires a recent version of the GEOmetadb SQLite file.  Download the file, and unzip it in the [data/](data/) directory:  
```R
library(GEOmetadb)
getSQLiteFile(destdir = "data/", destfile = "GEOmetadb.sqlite.gz")
```


### Run example metaGEO shiny-app:
```R
mypackage::runExample()
```

