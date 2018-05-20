library(GenomicDataCommons)
manifest <- read.table("G:/biodata/TCGA/gdc_manifest.2018-05-10.txt", header = TRUE, stringsAsFactors = FALSE)

# library(GenomicDataCommons)
# library(magrittr)
file_uuids <- manifest$id
head(file_uuids)

TCGAtranslateID = function(file_ids, legacy = TRUE) {
    info = files(legacy = legacy) %>%
        filter( ~ file_id %in% file_ids) %>%
        select('cases.samples.submitter_id') %>%
        results_all()
    # The mess of code below is to extract TCGA barcodes
    # id_list will contain a list (one item for each file_id)
    # of TCGA barcodes of the form 'TCGA-XX-YYYY-ZZZ'
    id_list = lapply(info$cases,function(a) {
        a[[1]][[1]][[1]]})
    # so we can later expand to a data.frame of the right size
    barcodes_per_file = sapply(id_list,length)
    # And build the data.frame
    return(data.frame(file_id = rep(ids(info),barcodes_per_file),
                      submitter_id = unlist(id_list)))
}

res = TCGAtranslateID(file_uuids)
head(res)
