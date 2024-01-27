# convert between symbol id and ensembl id

idconv <- function(from_id, from_type = "SYMBOL", to_type = "GENEID", db = "mouse", ...){
    require(EnsDb.Mmusculus.v79)
    require(EnsDb.Hsapiens.v86)
    require(AnnotationHub)
    require(magrittr)
    # ah <- AnnotationHub(localHub = TRUE)
    
    if(db == "mouse"){
        ensdb <- EnsDb.Mmusculus.v79
        # ensdb = ah[["AH109655"]]
    }else if(db == "human"){
        ensdb <- EnsDb.Hsapiens.v86
        # ensdb = ah[["AH109606"]]
    }
    
    
    mapIds(x = ensdb, 
           keys = from_id, 
           keytype = from_type, 
           column = to_type) 
}
