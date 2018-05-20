setClass("NeoSets",
         representation(
             summaryData = "ANY",
             geneExprs = "data.frame",
             geneCNV = "data.frame",
             geneMutation = "data.frame",
             sampleIdentifier = "character",
             sampleNeoSummary = "data.frame",
             sampleTIL = "data.frame"
         ),
         contains = c("MAF","GISTIC"), prototype = prototype(tidyData=NULL,
                                                 sampleIdentifier="Tumor_Sample_Barcode")
         )
## Set validity for NeoSets objects
# setValidity("NeoSets",
#             function(object){
#                 !is.null(object@sampleIdentifier)
#             })

testSets <- new("NeoSets",
                geneCNV = laml.plus.gistic@data,
                geneMutation = laml.plus.gistic@data,
                sampleNeoSummary = laml.plus.gistic@data,
                data = laml.plus.gistic@data,
                variants.per.sample = laml.plus.gistic@variants.per.sample,
                variant.type.summary = laml.plus.gistic@variant.type.summary,
                variant.classification.summary = laml.plus.gistic@variant.classification.summary,
                gene.summary = laml.plus.gistic@gene.summary,
                summary = laml.plus.gistic@summary,
                maf.silent = laml.plus.gistic@maf.silent,
                clinical.data = data.table::data.table(laml@clinical.data),
                cnv.summary = laml.gistic@cnv.summary,
                cytoband.summary = laml.gistic@cytoband.summary,
                cnMatrix = laml.gistic@cnMatrix,
                numericMatrix = laml.gistic@numericMatrix,
                gis.scores = laml.gistic@gis.scores,
                classCode = laml.gistic@classCode)


