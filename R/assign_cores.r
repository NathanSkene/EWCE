#' Assign cores
#'
#' Assign cores automatically for parallel processing, while reserving some.
#'
#' @param worker_cores Number (>1) or proportion (<1) of worker cores to use.
#' @param verbose Print messages.
#'
#' @return List of core allocations.
#'
#' @importFrom DelayedArray setAutoBPPARAM
#' @importFrom parallel detectCores
#' @keywords internal
assign_cores <- function(worker_cores = .90,
                         verbose = TRUE) {
    # Enable parallelization of HDF5 functions
    ## Allocate ~10% of your available cores to non-parallelized processes
    worker_cores <- if (is.null(worker_cores)) .90 else worker_cores
    total_cores <- parallel::detectCores()
    if (worker_cores < 1) {
        reserved_cores <- ceiling(total_cores * (1 - worker_cores))
        workers <- total_cores - reserved_cores
    } else {
        workers <- worker_cores
        reserved_cores <- total_cores - workers
    }
    messager(
        "+ ", workers, " core(s) assigned as workers (",
        reserved_cores, " reserved)."
    )
    DelayedArray::setAutoBPPARAM(BiocParallel::MulticoreParam(workers))
    #### Not allowed to use internal functions ####
    # DelayedArray:::set_verbose_block_processing(verbose)
    return(list(
        worker_cores = workers,
        reserved_cores = reserved_cores,
        total_cores = total_cores
    ))
}
