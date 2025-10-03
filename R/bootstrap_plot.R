#' Create a multi-page PDF with each cell type on a separate page
#'
#' This function creates a PDF where each unique value in a specified faceting 
#' variable (typically "CellType") is plotted on a separate page. Each plot shows 
#' expression data with points colored by expression levels.
#'
#' @param gene_data A data frame containing gene expression data with columns for 
#'   "boot", "hit", and a faceting variable (typically "CellType").
#' @param facet_var Character string specifying the column name to use for 
#'   separating plots onto different pages (e.g., "CellType").
#' @param output_file Character string specifying the path to save the PDF file.
#' @param add_labels Character string specifying the column name to use for point labels,
#'   or NULL to disable labels (default: NULL).
#' @param verbose Logical indicating whether to print progress messages (default: TRUE).
#'
#' @return No return value, called for side effect of creating a PDF file.
#'
#' @details The function creates QQ plots of "boot" vs "hit" values, with points 
#'   colored by "hit" values using a reversed viridis color scale. Each plot includes
#'   a title showing the cell type being displayed. If add_labels is provided, text labels
#'   from the specified column will be added to the points using ggrepel.
#'
#' @examples
#' \dontrun{
#' # Create plots without labels
#' .create_multipage_plot(
#'   gene_data = my_gene_data,
#'   facet_var = "CellType",
#'   output_file = "cell_type_plots.pdf"
#' )
#'
#' # Create plots with gene symbol labels
#' .create_multipage_plot(
#'   gene_data = my_gene_data,
#'   facet_var = "CellType",
#'   output_file = "cell_type_plots_labeled.pdf",
#'   add_labels = "gene_symbol"
#' )
#' }
#'
#' @import ggplot2
#' @importFrom viridis scale_color_viridis_c
#' @importFrom ggrepel geom_text_repel
#'
#' @keywords internal 
.create_multipage_plot <- function(gene_data, facet_var, output_file, 
                                   add_labels = NULL, verbose = TRUE) {
  # Check for required packages
  if(!is.null(add_labels) && !requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is needed for adding labels. Please install it.",
         call. = FALSE)
  }
  
  if (verbose) messager("Saving plot -->", output_file)
  # Check if the specified label column exists
  if(!is.null(add_labels) && !(add_labels %in% colnames(gene_data))) {
    stop(paste0("Label column '", add_labels, "' not found in data."),
         call. = FALSE)
  }
  
  # Get unique cell types
  cell_types <- unique(gene_data[[facet_var]])
  
  # Start PDF device
  pdf(output_file, width = 5, height = 5)
  
  # Create and print a plot for each cell type
  for (cell_type in cell_types) {
    # Subset data for this cell type
    subset_data <- gene_data[gene_data[[facet_var]] == cell_type, ]
    
    # Create base plot for this cell type
    g <- ggplot(subset_data, aes(x = boot, y = hit, color = hit)) +
      geom_point(size = 1, alpha = .75) + 
      xlab("Mean Bootstrap Expression") +
      ylab("Expression in cell type (%)\n") +
      scale_color_viridis_c(direction = -1) +
      ggtitle(paste("Cell Type:", cell_type)) +
      geom_abline(
        intercept = 0, 
        slope = 1, 
        linetype = "dashed",
        colour = ggplot2::alpha("red",.5)
      ) +
      theme_classic()
    
    # Add labels if a column name was provided
    if(!is.null(add_labels)) {
      g <- g + ggrepel::geom_text_repel(
        aes(label = .data[[add_labels]]),  # Use the column name dynamically
        size = 3,
        force_pull = 0.1, 
        box.padding = 0.35,
        point.padding = 0.5,
        segment.color = "grey50",
        max.overlaps = 25
      )
      
      if (verbose) messager("Added labels from column:", add_labels)
    }
    
    # Print the plot (adds it to the PDF)
    print(g)
    
    if (verbose) messager("Added plot for", cell_type)
  }
  
  # Close the PDF device
  dev.off()
  
}


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
    gene_data$symLab <- ifelse(gene_data$hit - gene_data$boot > 1, gene_data$gene, "") 
    
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
    #### Plot 1: Plot without gene names  ####
    .create_multipage_plot(
      gene_data = gene_data,
      facet_var = "CellType",  # Column name containing cell types
      output_file = files[[1]],
      verbose = verbose
    )
    
    #### Plot 2: Plot with gene names  ####
    .create_multipage_plot(
      gene_data = gene_data,
      facet_var = "CellType",  # Column name containing cell types
      output_file = files[[2]],
      add_labels = "symLab",
      verbose = verbose
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
        # cat("Names of melt boot:", paste(names(melt_boot), collapse = " "))
        # wd <- 1 + length(unique(melt_boot[, hit])) * 0.25 
        plots[["plot3"]] <- g3
        messager("Saving plot -->", files[[3]], v=verbose)
        ggplot2::ggsave(filename = files[[3]], 
                        plot = g3,
                        width = 8, 
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
        # wd <- 1 + length(unique(melt_boot[,4])) * 0.25 
        plots[["plot4"]] <- g4
        messager("Saving plot -->", files[[4]], v=verbose)
        ggplot2::ggsave(filename = files[[4]], 
                        plot = g4,
                        width = 8,
                        height = 4) 
    }
   
    if(isTRUE(show_plot)) methods::show(plots)
    return(list(plots=plots,
                paths=files))
}


