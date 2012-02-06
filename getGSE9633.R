## THIS IS MY NEW FILE
#####

options(stringsAsFactors=F)

require(synapseClient)
require(metaGEO)

## LOOK FOR GSE9633
myQuery <- synapseQuery("SELECT * FROM dataset WHERE name == 'GSE9633'")
gse9633 <- getEntity(myQuery$dataset.id[myQuery$dataset.parentId==102611])
gse9633layers <- getDatasetLayers(gse9633)
expr9633 <- loadEntity(gse9633layers$layer.id)$objects$data ## metaGEO class

exprMat <- expr9633@assayData[["exprs"]]

