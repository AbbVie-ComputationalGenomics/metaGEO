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


### Run example metaGEO shiny-app
Copy or move the `.Rdata` file generated from the vignette to the shiny app's data location:
```bash
cp vignettes/ObjectWithGeoData_psoriasis.Rdata inst/shiny-examples/data/
cd inst/shiny-examples/

Start the app
```R
library(shiny)
runApp("metaGEO")
```

