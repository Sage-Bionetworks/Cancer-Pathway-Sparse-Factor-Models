## THIS IS MY NEW FILE
#####

options(stringsAsFactors=F)

require(synapseClient)
require(metaGEO)

## LOOK FOR GSE5851
myQuery <- synapseQuery("SELECT * FROM dataset WHERE name == 'GSE5851'")
gse5851 <- getEntity(myQuery$dataset.id[myQuery$dataset.parentId==102611])
gse5851layers <- getDatasetLayers(gse5851)
expr5851 <- loadEntity(gse5851layers$layer.id)$objects$data ## metaGEO class

exprMat <- expr5851@assayData[["exprs"]]

