#' allocation
#'
#' Allocates changes in land use demand across the landscape
#'
#' @param lu Matrix containing fractional land use data. Columns are classes, rows are cells. Each row sums to 1.
#' @param sm A matrix containing the predicted land use suitability model for the current time step.
#' @param params A list of model parameters _resolution_, _max_dev_, _growth_, _no_change_ (see details).
#' @param dmd Matrix containing demand changes to be allocated. The first row is the current land use supply, the second row the land use demand for the next time step. Output of `fliutes::demand`.
#' @param ln A matrix containing neighbourhood variables of the current time step. These can be estimated by using the `flutes::neighbourhood` function.
#' @param constraint `TRUE` or `FALSE`. When true, constraints are applied (see details).
#' @param pa default is `pa = NULL` or raster layer containing protected areas (where the raster contains 0) from which changes are masked. All other areas (where the raster is 1) are allowed to change.
#' @return A matrix containing predicted land use fractions (same format as input matrix)
#'
#' @details Model parameters that can be specified with `params` argument: _max_dev_ is the % deviation allowed between allocated supply and prescribed demand. _resolution_ is the integer count to which fractions are converted for multinomial model. i.e. a fraction of 0.2 would become round(0.2 * resolution). _growth_ is the % of the landscape (in terms of cells) where land use was 0 in a land use and contained a fraction of that land-use in the next time step, averaged across the observed time series. _no_change_ is a vector containing the indices of land use types that are not allowed to change at all. Applied constraints when `constraint = TRUE`, most cells on which land use in a class is 0 are masked from increases in that class to the next time step. Only a small subset of cells for which this is the case are allowed to change. The number of these cells is determined by _growth_, which can be estimated from historic data, literature review or expert advice.
#' @import wrswoR
#' @export


allocation <- function(lu, sm, params, dmd, ln, constraint, pa = NULL){

  #number of land use classes and number of cells
  K <- ncol(lu)
  n <- nrow(lu)

  resolution <- params$resolution
  max_dev <- params$max_dev
  max_iter <- params$max_iter
  growth <- params$growth
  no_change <-  1:K%in%params$no_change

  #Turn intiital land use into integer counts
  #p_t0 <- integerify(lu, resolution = resolution)
  p_t0 <- round(lu * resolution)

  #Convert demand time series to integer counts
  supply_t0 <- colSums(p_t0) #Initial land use
  demand_traj <- integerify(x = dmd, resolution = sum(supply_t0))
  demand_t1 <- demand_traj[2,] #Second row of demand trajectory is demand to be allocated
  supply_t1_candidate <- demand_traj[1,] #First time step land use supply becomes "candidate" that gets recalculated in each iteration until it meets demand_t1

  dev_diff <- abs(diff(demand_traj))/demand_traj[1,] * 100 #Initital candidate % deviation of current supply from demand, gets recalculated in each iteration until it's below max_dev
  dev_diff[which(is.na(dev_diff))] <- 0

  #Counter
  count <- 0

  if(!is.null(pa)){
    pa_inds <- which(pa==0)
    supply_t0_pa <- colSums(p_t0[pa_inds,])
    nopa_inds <- which(pa==1)
    p_t0 <- p_t0[nopa_inds,]
    ln <- ln[nopa_inds,]
    sm <- sm[nopa_inds,]
  }

  #First land use map becomes "candidate" on which allocations take place iteratively
  p_t1_candidate <- p_t0

  #Determine which classes are 0 on each cell (so we can constrain growth to those that aren't 0)
  if(constraint){
    are_zero <- p_t0 == 0
    inds_list <- list()
    i <- 1
    for(i in 1:K){
      inds_zero_all <- which(are_zero[,i]) #all that are 0 should stay 0, except for...
      inds_zero <- inds_zero_all[which(ln[inds_zero_all,i]!=0)] # a subset of cells near where land use already exists (neigbourhood)
      size <- length(are_zero[,i]) * ((growth[i])/100) # we want to sample this many from that subset to remain.
      if(size < length(inds_zero)){
        keep <- inds_zero[sample_int_rank(length(inds_zero), size = size, prob = sm[inds_zero,i])]
        inds_list[[i]] <- inds_zero_all[-which(inds_zero_all%in%keep)]
      }
      if(size > length(inds_zero)){
        leftover <- size - length(inds_zero)
        inds_leftover <- which(!inds_zero_all%in%inds_zero)
        keep <- inds_leftover[-sample_int_rank(length(inds_leftover), size = min(c(leftover, length(inds_leftover))), prob = sm[inds_leftover,i])]
        inds_list[[i]] <- inds_zero_all[-which(inds_zero_all%in%c(inds_zero, keep))]
      }
    }
  }

  # Iterative allocation
  while (any(dev_diff > max_dev) & count < max_iter) {

    #Counter increment
    count <- count + 1

    #Demand to be allocated
    demand_change <- demand_t1 - supply_t1_candidate

    if(!is.null(pa)){
      demand_change <- demand_change - supply_t0_pa
    }

    #Calculate change factor (a),by how much do we have to multiply the candidate land use proportions to satisfy the modelled suitability.
    ideal_change <- (sm * resolution) / p_t1_candidate

    both_0 <- is.na(ideal_change) #this determines which cells are 0 in the suitability map and the current iteration. We keep those as they are.
    ideal_change[both_0] <- 1

    cand_0 <- !is.finite(ideal_change) #This determines which cells are 0 in p_t1_candidate
    ideal_change[cand_0] <- sm[cand_0]

    #Calculate Relative suitability (r) from change factors. Sums to 1 in each column
    rel_suitability <- ideal_change %*% diag(1/colSums(ideal_change))

    #Allocate demand change between pixels (d), wieghted by r (d * r = m) and adjusted by stepsize (default is 1 but can be smaller for more fine-scale allocations: more stable, but slower)
    target_lu_change_pixel <-  rel_suitability %*% diag(demand_change)

    #Add changes to candidate map and make everyting positive.
    p_t1_proposal <- p_t1_candidate + target_lu_change_pixel
    p_t1_proposal <- pmax(p_t1_proposal, 0)

    #Where land use demand is 0, make all cells 0 in that land use
    if(any(demand_t1 == 0)){
      p_t1_proposal[, which(demand_t1 == 0)] <- 0
    }

    if(constraint){
      for(i in (1:K)){
        p_t1_proposal[inds_list[[i]], i] <- 0
      }
    }

    #Turn proposed land use map into integers
    p_t1_candidate <- integerify(x = p_t1_proposal,  resolution = resolution, no_decrease = no_change, z = p_t1_candidate)

    #ii Calcuate new candidate supply, i.e. the supply of the currently proposed candidate
    supply_t1_candidate <- colSums(p_t1_candidate)

    #Recalculate % deviation of candidate supply from demand
    diff <- abs(demand_t1 - supply_t1_candidate)
    if(!is.null(pa)){
      diff <- abs(demand_t1 - (supply_t1_candidate + supply_t0_pa))
    }
    dev_diff <- diff/demand_t1 * 100
    dev_diff[which(is.na(dev_diff))] <- 0

    cat("\r", paste0("Iteration: ", count, "    "), "Deviation from target per class [%]: ", paste(round(dev_diff, 3), sep = " "))
  }

  #When allocations are ready, return result
  if(!is.null(pa)){
    pred_out <- matrix(NA, ncol = K, nrow = n)
    pred_out[nopa_inds,] <- p_t1_candidate/resolution
    pred_out[pa_inds,] <- lu[pa_inds,]
  }else{
    pred_out <- p_t1_candidate/resolution
  }
  colnames(pred_out) <- colnames(lu)
  pred_out
}
