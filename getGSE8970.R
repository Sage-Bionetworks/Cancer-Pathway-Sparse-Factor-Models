## THIS IS MY NEW FILE
#####

options(stringsAsFactors=F)

require(synapseClient)
require(metaGEO)

## LOOK FOR GSE8970
myQuery <- synapseQuery("SELECT * FROM dataset WHERE name == 'GSE8970'")
gse8970 <- getEntity(myQuery$dataset.id[myQuery$dataset.parentId==102611])
gse8970layers <- getDatasetLayers(gse8970)
expr8970 <- loadEntity(gse8970layers$layer.id)$objects$data ## metaGEO class

exprMat <- expr8970@assayData[["exprs"]]

