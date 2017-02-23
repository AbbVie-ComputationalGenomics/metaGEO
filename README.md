# CompGenomics-metaGEO

### Clone this repository
```bash
git clone git@pig.abbvienet.com:grundaj/CompGenomics-metaGEO.git
```
### Start an R session

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
          'devtools'))
```

## Run the *psoriasis* vignette

```R
$ setwd(CompGenomics-metaGEO/vignettes)
```

Perform the initial rendering:
```R
library(rmarkdown)
render('build_metageo_instance.Rmd')
```

Upon successful rendering, load the resulting `build_metageo_instance.html` file in your browser.


### Run example metaGEO shiny-app
Copy or move the `.Rdata` file generated from the vignette to the shiny app's data location:
```bash
$ cp vignettes/ObjectWithGeoData_psoriasis.Rdata inst/shiny-examples/data/
$ cd inst/shiny-examples/
```

Start the app
```R
> library(shiny)
> runApp("metaGEO")
```

