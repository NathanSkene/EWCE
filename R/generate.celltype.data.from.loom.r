generate_celltype_data_from_loom <- function(loomFilePath,lvl1_label="level1class",lvl2_label="level2class",no_cores=NA){
  cell.attrs = get_loom_cell_attrs_as_df(loomFilePath)
  
  lfile <- loomR::connect(filename = loomFilePath, mode = "r+")
  #get_ct_means_lvl1 = function(x){
  #  ct_means = t(aggregate(x,by=list(cell.attrs[[lvl1_label]]),FUN=mean)[,-1])
  #  return(t(ct_means))
  #}
  get_ct_means_lvl2 = function(x){
    ct_means = t(aggregate(x,by=list(cell.attrs[[lvl2_label]]),FUN=mean)[,-1])
    return(t(ct_means))
  }

  map_Mean_in_CT_lvl2 <- lfile$map(FUN = get_ct_means_lvl2, MARGIN = 1, chunk.size = 500, dataset.use = "matrix", display.progress = TRUE) #,index.use=1:600) # lengt9430\
  rownames(map_Mean_in_CT_lvl2) = sort(unique(cell.attrs[[lvl2_label]]))
  colnames(map_Mean_in_CT_lvl2) = lfile$row.attrs[["HGNC"]][] #gene.attrs[["HGNC"]]
  save(map_Mean_in_CT_lvl2,file=sprintf("Data/Tidy/map_Mean_in_CT_%s.rda",lvl2_label))
  
  
  lvls = unique(cell.attrs %>% select(lvl1_label,lvl2_label))
  mm_t = data.frame(map_Mean_in_CT_lvl2)
  mm_t[,lvl2_label] = rownames(mm_t)
  mm_t2 = merge(mm_t,lvls,by=lvl2_label)
  mm_t3 = mm_t2 %>% select(-c(lvl2_label))
  #lvl1_means = plyr::ddply(mm_t3,.(level1class_dx),colwise(mean))
  
  
  mm_t3[,"sub_level"] = as.character(mm_t3[,lvl1_label])
  mm_t3 = mm_t3 %>% select(-c(lvl1_label))
  
  
  #mm_t3[,lvl1_label] = gsub(" ","____",mm_t3[,lvl1_label])
  
  library(dplyr)
  library(plyr)
  lvl1_means = plyr::ddply(mm_t3,.(sub_level),colwise(mean))
  
  #lvl1_means = plyr::ddply(mm_t3,mm_t3$sub_level,colwise(mean))
  #lvl1_means = plyr::ddply(mm_t3 %>% select(-level1class),as.character(mm_t3[,lvl1_label]),colwise(mean))
  
  
  #lvl1_means = plyr::ddply(mm_t3,mm_t3[,lvl1_label],colwise(mean))
  # lvl1_means$level1class_dx = gsub("\\^\\^\\^\\^"," ",lvl1_means$level1class_dx)
  rownames(lvl1_means) = lvl1_means$sub_level
  lvl1_means = lvl1_means %>% select(-sub_level)
  map_Mean_in_CT_lvl1 = data.matrix(lvl1_means)
  colnames(map_Mean_in_CT_lvl1) = colnames(map_Mean_in_CT_lvl2)
  
  #map_Mean_in_CT_lvl1 <- lfile$map(FUN = get_ct_means_lvl1, MARGIN = 1, chunk.size = 500, dataset.use = "matrix", display.progress = TRUE) #,index.use=1:600) # lengt9430\
  #rownames(map_Mean_in_CT_lvl1) = sort(unique(cell.attrs[[lvl1_label]]))
  #colnames(map_Mean_in_CT_lvl1) = lfile$row.attrs[["HGNC"]][] #gene.attrs[["HGNC"]]
  save(map_Mean_in_CT_lvl1,file=sprintf("Data/Tidy/map_Mean_in_CT_%s.rda",lvl1_label))
  
  lfile$close_all()
  
  calculate.specificity.for.level <- function(ctd_oneLevel){
    normalised_meanExp = t(t(ctd_oneLevel$mean_exp)*(1/colSums(ctd_oneLevel$mean_exp)))
    ctd_oneLevel$specificity = normalised_meanExp/(apply(normalised_meanExp,1,sum)+0.000000000001)
    return(ctd_oneLevel)
  }
  
  ctd = list()
  ctd[[1]] = list()
  ctd[[1]]$mean_exp = t(map_Mean_in_CT_lvl1)
  ctd[[2]] = list()
  ctd[[2]]$mean_exp = t(map_Mean_in_CT_lvl2)  
  
  if(is.na(no_cores)){no_cores = detectCores()-1}
  cl <- makeCluster(no_cores,type="FORK")  
  ctd = parallel::mclapply(ctd,calculate.specificity.for.level,mc.cores=no_cores,mc.preschedule = FALSE)
  
  crapGenes = rowSums(ctd[[2]]$mean_exp<0.02 & ctd[[2]]$specificity>0.04)
  keepGenes = setdiff(rownames(ctd[[2]]$mean_exp),names(crapGenes)[crapGenes>=1])
  # ctd[[1]]$mean_exp = ctd[[1]]$mean_exp[gsub("-","\\.",keepGenes),]
  # ctd[[2]]$mean_exp = ctd[[2]]$mean_exp[gsub("-","\\.",keepGenes),]
  # ctd[[1]]$specificity = ctd[[1]]$specificity[gsub("-","\\.",keepGenes),]
  # ctd[[2]]$specificity = ctd[[2]]$specificity[gsub("-","\\.",keepGenes),]
  ctd[[1]]$mean_exp = ctd[[1]]$mean_exp[keepGenes,]
  ctd[[2]]$mean_exp = ctd[[2]]$mean_exp[keepGenes,]
  ctd[[1]]$specificity = ctd[[1]]$specificity[keepGenes,]
  ctd[[2]]$specificity = ctd[[2]]$specificity[keepGenes,]
  
  colnames(ctd[[1]]$mean_exp) = trimws(colnames(ctd[[1]]$mean_exp))
  colnames(ctd[[1]]$specificity) = trimws(colnames(ctd[[1]]$specificity))
  
  
  # ADD DENDROGRAM DATA TO CTD
  ctd = lapply(ctd,EWCE::bin.specificity.into.quantiles,numberOfBins=40)
  library(ggdendro)
  library(EWCE)
  ctd = lapply(ctd,prep.dendro)
  return(ctd)
}

