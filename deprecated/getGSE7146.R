## THIS IS MY NEW FILE
#####

options(stringsAsFactors=F)

require(synapseClient)
require(metaGEO)

## LOOK FOR GSE7146
myQuery <- synapseQuery("SELECT * FROM dataset WHERE name == 'GSE7146'")
gse7146 <- getEntity(myQuery$dataset.id[myQuery$dataset.parentId==102611])
gse7146layers <- getDatasetLayers(gse7146)
expr7146 <- loadEntity(gse7146layers$layer.id[gse7146layers$layer.platform=="hgu133a"])$objects$data ## metaGEO class

exprMat <- expr7146@assayData[["exprs"]]

