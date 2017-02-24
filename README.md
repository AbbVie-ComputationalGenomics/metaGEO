# CompGenomics-metaGEO
## Quickstart guide
Following this guide will build an example metaGEO app comparing multiple psoriasis studies.  More details on how to create a new metaGEO can be found by following the vignette.
### Clone this repository
```bash
git clone git@pig.abbvienet.com:grundaj/CompGenomics-metaGEO.git
```
### Start an Rstudio session

### R package requirements
To install all dependencies run:
```  
source("https://bioconductor.org/biocLite.R")
  biocLite(c('shiny',
          'shinyjs',
          'shinyBS',
          'lattice',
          'logging',
          'hgu133plus2.db',
          'org.Hs.eg.db', 
          'knitr',
          'gdata',
          'dbConnect',
          'RSQLite',
          'GEOmetadb',
          'xtable',
          'rmarkdown',
          'devtools',
          'reshape2'))
```

### Pandoc requirement
Please install the appropriate system build of Pandoc

[https://github.com/jgm/pandoc/releases/tag/1.19.2.1]

## Run the *psoriasis* vignette

```R
setwd('CompGenomics-metaGEO/vignettes')
```

Perform the initial rendering (May take a few minutes depending on network speed):
```R
library(rmarkdown)
render('build_metageo_instance.Rmd')
```

Upon successful rendering, load the resulting `build_metageo_instance.html` file in your browser.


### Run example metaGEO shiny-app
Copy or move the `.Rdata` file generated from the vignette to the shiny app's data location:
```R
file.rename('ObjectWithGeoData_psoriasis.Rdata', '../inst/shiny-examples/data/ObjectWithGeoData_psoriasis.Rdata')
setwd('../inst/shiny-examples/')
```

Start the app
```R
library(shiny)
runApp("metaGEO")
```

