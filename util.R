populateSCANExpressionSets<-function(in_SCAN, in_RMA){
  
  ############################################################################
  # this step is now filling in all of the proper slots of the ExpressionSet
  # object created by SCAN.UPC
  # we do this by borrowing slots that already have been created from the
  # raw GeneFeatureSet object and the RMA normalized ExpessionSet object
  ############################################################################
  
  # check that the order of the column names in normalized is the same as in affyGeneFS before mapping
  #all(pData(affyGeneFS)$CELname==colnames(normalized))
  if(all(pData(in_RMA)$CELname==colnames(in_SCAN))==FALSE) 
    stop("Order of column names of SCAN normalized data does not match order of column names of RMA normalized data. PLease check your data")
  #must be TRUE
  
  phenoData(in_SCAN) <- phenoData(in_RMA)
  varMetadata(in_SCAN)<-varMetadata(in_RMA)
  
  #phenoData(normalized) <- phenoData(affyGeneFS)
  #varMetadata(normalized)<-varMetadata(affyGeneFS)
  
  #all(colnames(normalized)==row.names(pData(normalized)))
  if(all(colnames(in_SCAN)==row.names(pData(in_SCAN)))==FALSE)
    stop("Check your data because the row.name order in pData of SCAN normed data doesn't match column name order in SCAN normed data")
  
  colnames(exprs(in_SCAN))<-row.names(pData(in_SCAN))
  
  #dimnames(se.exprs(normalized)) = dimnames(exprs(normalized))
  
  #add featuredata, must come from geneCoreRMA because affyGeneFS doesn't have the same level of features
  
  featureData(in_SCAN)  <- featureData(in_RMA)
  protocolData(in_SCAN) <- protocolData(in_RMA)
  annotation(in_SCAN)   <- annotation(in_RMA)
  # pData(in_SCAN)<-cbind(pData(in_SCAN),protocolData(in_SCAN)@data)
  
  return(in_SCAN)
}

#populateSCANExpressionSets(normalized, geneCoreRMA)