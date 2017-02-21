# CompGenomics-metaGEO

### R package requirements
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

## Run the *psoriasis* vignette
```bash
$ cd vignettes
$ ls
build_metageo_instance.Rmd  psoriasis_info.R
```

Start an R session, and perform the initial rendering:
```R
> library(rmarkdown)
> rmarkdown::render('build_metageo_instance.Rmd')
```

Upon successful rendering, load the resulting `build_metageo_instance.html` file in your browser.


### Run example metaGEO shiny-app: (in progress)
```R
mypackage::runExample()
```

