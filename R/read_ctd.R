# read_ctd <- function(CTD_meta_row){
#     sourceRData <- function(fileName){
#         repmis::source_data(fileName)
#         get(ls()[ls() != "fileName"])
#     }
#     ROW <- CTD_meta_row
#     if(endsWith(tolower(ROW$url),".rds")){
#         print("Reading in RDS file...")
#         ctd <- readRDS(url(ROW$url,"rb"))
#     }else {
#         print("Reading in RDA file...")
#         ctd <- sourceRData(ROW$url)
#     }
#     return(ctd)
# }
