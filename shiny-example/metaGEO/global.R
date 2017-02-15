packagedir <- '../'
print('loading Rdata')
load(paste0(packagedir, 'data/ObjectWithGEOdata_psoriasis.Rdata'))
metageo_title <- "Psoriasis"

jsCode <- "
shinyjs.flip_checkbox = function (params) {
  // cb_value and onoff should be passed as named arguments
  var defaultParams = {
    cb_value: 'none',
    onoff: 'off'
  }
  params = shinyjs.getParams(params, defaultParams);
  cb_value = params.cb_value;
  onoff = params.onoff;
  console.log('called flip_checkbox for: ' + cb_value + '  state: ' + onoff);
  var cbs = document.querySelectorAll(\"input[value=\"+cb_value+\"]\"); 
  for(var i = 0; i < cbs.length; i++) {
    if(onoff == 'off') {
      cbs[i].checked = false;
    }else{
      cbs[i].checked = true;
    }
  }
}

shinyjs.click_generatePlots = function () {
  document.getElementById(\"generatePlots\").click();
}
"

genes_list <- read.table(paste0(packagedir, 'data/hgnc_complete_set.txt.gz'), as.is=T, sep='\t', 
                         quote='', comment.char = '', header=T)[,2]

datasets <- names(data)

AtoB <- c()
geneindex <- ''

plotmethods <- c('Dot', 'Box', 'Violin')
textmethods <- c('Long', 'Ragged', paste('Grouping', datasets, sep='_'))
