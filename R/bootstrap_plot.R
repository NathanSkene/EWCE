#' Bootstrap plot
#' 
#' Plot bootstrap enrichment results. 
#' Support function for \link[EWCE]{generate_bootstrap_plots}. 
#' @param gene_data Output from \link[EWCE]{compute_gene_scores}. 
#' @param save_dir Directory to save plots to.
#' @param listFileName listFileName
#' @param exp_mats Output of \code{generate_bootstrap_plots_exp_mats}.
#' @param show_plot Print the plot.
#' @param signif_ct Significant celltypes to include the plots.
#' @inheritParams ggplot2::facet_grid
#' @returns Null output.
#' 
#' @keywords internal 
#' @importFrom data.table .I
bootstrap_plot <- function(gene_data,
                           exp_mats=NULL,  
                           save_dir = file.path(tempdir(),"BootstrapPlots"),
                           listFileName,
                           signif_ct=NULL,
                           hit_thresh = 25,
                           facets = "CellType",
                           scales = "free_x",
                           show_plot = TRUE,
                           verbose = TRUE) {
    # devoptera::args2vars(bootstrap_plot)
    
    requireNamespace("ggplot2")
    requireNamespace("patchwork")
    requireNamespace("ggrepel")
    Pos <- Rep <- Exp <- p <- significant <- CellType <- NULL;
    
    exp_mats_msg <- paste(
        "Cannot create bootstrap distribution plots",
        "without exp_mats"
    )
    plots <- list()
    #### Set up save path ####
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)  
    gene_data$symLab <- ifelse(gene_data$hit > 25, gene_data$gene, "") 
    
    #### Filter gene_data ####
    if(!is.null(signif_ct)){
        gene_data <- gene_data[CellType %in% signif_ct,]
    }
    #### Prepare file paths ####
    files <- file.path(save_dir,
                       sprintf(c("qqplot_noText____%s.pdf",
                                 "qqplot_wtgene____%s.pdf",
                                 "bootDists____%s.pdf",
                                 "bootDists_LOG____%s.pdf"),
                               listFileName)
              )
    for(f in files){
        dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)   
    }
    ## Plot several variants of the graph ##
    add_line <- function(){
      geom_abline(
        intercept = 0, 
        slope = 1, 
        linetype = "dashed",
        color = "red", 
        alpha = 0.5
      ) 
    } 
    #### Plot 1: Plot without gene names  ####
    g1 <- ggplot(gene_data, aes(x = boot, y = hit,  color = hit)) +
        geom_point(size = 1,  alpha = .75) + 
        xlab("Mean Bootstrap Expression") +
        ylab("Expression in cell type (%)\n") +
        scale_color_viridis_c() +
        facet_grid(facets = facets, scales = scales) + 
        add_line() +
        theme_classic() 
    plots[["plot1"]] <- g1 
    messager("Saving plot -->", files[[1]], v=verbose)
    ggsave(
      filename = files[[1]], 
      plot = g1,
      width = 4, 
      height = 3.5
    )
    
    #### Plot 2: Plot with gene names  ####
    g2 <- g1 + 
        geom_text_repel(
          mapping = aes(label = symLab), 
          alpha = 0.75,
          segment.alpha = 0.75,
          max.overlaps = 25,
          force_pull = 0.5
        ) # +
        # scale_x_discrete(expand = expansion(mult = c(0,.15))) +
        # scale_y_discrete(expand = expansion(mult = c(0,.15))) 
    plots[["plot2"]] <- g2 
    messager("Saving plot -->", files[[2]], v=verbose)
    ggsave(
      filename = files[[2]], 
      plot = g2,
      width = 4, 
      height = 3.5
    ) 
   
 
    #### Plot 3 ####
    if(is.null(exp_mats)){
        messager(exp_mats_msg)
        files <- files[seq_len(2)]
    } else { 
        melt_boot <- 
            (
                lapply(exp_mats, function(x){
                    data.table::as.data.table(x)[,Rep:=paste0("rep",.I)]
                }) |>
                    data.table::rbindlist(use.names = TRUE, 
                                          idcol = "CellType") |>
                    data.table::melt(id=c("CellType","Rep"),
                                     variable.name="Pos",
                                     value.name = "Exp") 
            )[,Pos:=factor(as.integer(gsub("V","",Pos)), 
                              ordered=TRUE)]
        levels(melt_boot$Pos) <- rev(levels(melt_boot$Pos) )
        melt_boot[,Exp:=Exp*100]
        #### Filter celltypes ####
        if(!is.null(signif_ct)){
            melt_boot <- melt_boot[CellType %in% signif_ct,]
        }
        # actVals_ordered <- data.table::merge.data.table(gene_order, actVals)
        # actVals_ordered[,pos:=factor(pos,levels = unique(sort(pos)), ordered = TRUE)]
        # Plot with bootstrap distribution
        g3 <- ggplot(melt_boot) +
            geom_boxplot(aes_string(x = "Pos", y = "Exp"),
                         outlier.size = 0) +
            geom_point(aes_string(x = "rank", y = "hit"),
                       color = ggplot2::alpha("red",.75), 
                       data = gene_data
            ) + 
            ylab("Expression in cell type (%)\n") +
            xlab(expression("Least specific" %->% "Most specific")) +
            scale_x_discrete(breaks = NULL) +
            theme_graph() +
            facet_wrap(facets=facets)
        print(sprintf("Names of melt boot: %s", paste(names(melt_boot), collapse = " ")))
        wd <- 1 + length(unique(melt_boot[, hit])) * 0.25 
        plots[["plot3"]] <- g3
        messager("Saving plot -->", files[[3]], v=verbose)
        ggplot2::ggsave(filename = files[[3]], 
                        plot = g3,
                        width = wd, 
                        height = 3.5) 
        
        
        #### Plot 4: Plot with bootstrap distribution ####
        melt_boot2 <- merge(melt_boot, gene_data, 
                            by.x = c("CellType","Pos"),
                            by.y = c("CellType","rank"))
        data.table::setkeyv(melt_boot2,c("CellType","hit")) 
        melt_boot2[,significant:=factor(p<0.05, levels = c(TRUE,FALSE), 
                                        ordered = TRUE)] 
        #### Filter celltypes ####
        if(!is.null(signif_ct)){
            melt_boot2 <- melt_boot2[CellType %in% signif_ct,]
        }
        plts <- lapply(unique(melt_boot2$CellType), function(cc){
            dat <- melt_boot2[CellType==cc] 
            dat$gene <- factor(dat$gene, unique(dat$gene), ordered=TRUE)
            ggplot(dat) +
                geom_boxplot(aes_string(x = "gene", y = "Exp", fill="significant", 
                                        color="significant"),
                             outlier.size = 0) +
                scale_color_manual(values=c(ggplot2::alpha("blue",1),
                                            "grey"), drop=FALSE) +
                theme_graph() +
                theme(axis.text.x = element_text(angle = 54, hjust = 1)) +
                geom_point(aes_string(x = "gene", y = "hit"),
                           color=ggplot2::alpha("red",.5)
                ) +
                scale_fill_manual(values=c(ggplot2::alpha("blue",.25),
                                           ggplot2::alpha("white",.25)
                                           ), 
                                  drop=FALSE) +
                ylab(NULL) +
                xlab(NULL) +
                ylim(c(0,100)) +
                facet_wrap(facets = facets, 
                           scales = scales)
        })
        g4 <- patchwork::wrap_plots(plts) + 
            patchwork:: plot_layout(guides = 'collect') +
            labs(x=expression("Least specific" %->% "Most specific"),
                 y="Expression in cell type (%)\n") 
        #### Save ####
        wd <- 1 + length(unique(melt_boot[,4])) * 0.25 
        plots[["plot4"]] <- g4
        messager("Saving plot -->", files[[4]], v=verbose)
        ggplot2::ggsave(filename = files[[4]], 
                        plot = g4,
                        width = wd, 
                        height = 4) 
    }
   
    if(isTRUE(show_plot)) methods::show(plots)
    return(list(plots=plots,
                paths=files))
}
