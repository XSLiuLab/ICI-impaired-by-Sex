# Source code from biostar and mofied for easier use
TCGAtranslateID = function(file_ids, legacy = TRUE) {
    # file_ids: id of download files in TCGA manifest file
    # legacy: if TRUE, use hg19
    if(!(require(GenomicDataCommons) & require(magrittr))){
        install.packages(c("GenomicDataCommons", "magrittr"), dependencies = TRUE)
    }
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


# manifest <- read.table("gdc_manifest_filePath", header = TRUE, stringsAsFactors = FALSE)
# file_uuids <- manifest$id
# head(file_uuids)
# res = TCGAtranslateID(file_uuids)
# head(res)
