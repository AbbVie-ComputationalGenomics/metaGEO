# data required for building metageo instance's Rdata object

data_filename <- 'ObjectWithGeoData_psoriasis.Rdata'

datasets <- c(
  'GSE13355',
  'GSE14905',
  'GSE57225'
)


# data for making grid of groupings.  Run the vignette to generate the unique column values for
# each GSE.  Identifying the groupings are very much a manual process.

#  The main idea:
#
#          group_name_1   group_name_2    group_name_3    ...
#       -----------------------------------------------
#  GSE1 |     regex         regex            XXX
#  GSE2 |      XXX          regex            XXX
#  GSE3 |      XXX          regex           regex

# group_names are the headings to the grouping grid
#group_names <- c('1. PBMCs Crohns', '2. PBMCs Normal', '3. PBMCs UC','4. Biopsy Crohns', '5. Biopsy UC', 
#                 '6. Biopsy healthy', '7. Biopsy other', '8. IECs Crohns', '9. IECs Healthy' )
group_names <- c('1. Psoriasis non-lesion', '2. Psoriasis lesion', '3. Disease Other', '4. Healthy')

# designate the annotation column that holds the metadata that will describe the samples by the above group_names. 
# position in the list correlates to the position in 'datasets', can be a vector or single column name.
#columnList<-list(
#  'description', #1
#  'description', #2 
#  'source_name_ch1', #3
#  'source_name_ch1', #4
#  c('source_name_ch1'), #5
#  )
#  etc...

columnList <- list(
  'characteristics_ch1', # 1
  c('source_name_ch1', 'characteristics_ch1'),
  'characteristics_name_ch2.1'
)


# specify a regex for each position in the grouping grid.  each row (vector) corresponds to a GSE
# use 'XXX' when a GSE does not contain data associated with that grouping/header.  When muultiple columns
# contain metadata of interest (vector in the column List), those usually require an or "|" operator in the regex
#regexList <- list(
#  c("PBMCs from Crohn's Disease subject","PBMCs from Normal subject","PBMCs from Ulcerative Colitis subject", rep('XXX', 6)),     #1
#  c(rep('XXX', 3), "biopsy from crohns disease subject",'XXX', "biopsy from healthy subject", rep('XXX', 3)), #2
#  c(rep('XXX', 3), "colonoscopic biopsy from adult with affected colon with Crohn.*s disease","colonoscopic biopsy from adult with affected colon from adult with UC","colonoscopic biopsy from normal adult","colonoscopic biopsy from adult with affected colon from adult with IC", rep('XXX', 2)), #3
#  etc...
regexList <- list(
  c('^uninvolved skin from cases', '^involved skin from cases', 'XXX', 'controls'),
  c('Skin biopsy, non-lesional', 'Skin biopsy, lesional', 'XXX', 'Skin biopsy, normal'),
  c('non-involved skin', 'Lesional skin, psoriasis', 'Lesional skin, eczema', 'XXX')
)

