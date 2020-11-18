#' Plot community detections across site
#'
#' Visualize results of community detection output from
#' \link[FCharmony]{test_communities}
#'
#' @param x List of community detection results
#' @param names 
#' @param which
#' @param atlas
#' @param atlas.lab
#' @param match.labs Whether to modify community labels to match across levels.
#' @param overlap.obj If matching labels, which objective function to use. All
#'   are sum over all ROIs. 1 is number of overlapping labels, 2 is whether or
#'   not all labels match, 3 is the number of unique labels (minimize), 4 is
#'   a weighted sum between 1 and 2, giving more weight to having all same. 5 is
#'   a weighted sum between 1 and 3, penalizing more distinct labels.
#' @param manual.labs Optional vector of length equal to the number of
#'   communities. Recodes community numbers. For example, `c(3,2,1)` swaps the
#'   labels of communities 1 and 3.
#'
#' @return
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom RColorBrewer brewer.pal
#' @export
#'
#' @examples
plot_matlab_comm <- function(x,
                             grps,
                             which = 1,
                             atlas = NULL,
                             atlas.lab = NULL,
                             match.labs = FALSE,
                             overlap.obj = 1,
                             manual.labs = NULL,
                             colors = NULL) {
  comms <- t(x[[which]][[1]])
  colnames(comms) <- grps

  ## match community labels based on overlap
  switch(
    overlap.obj,
    "1" = {overlap_obj <- function(x) max(tabulate(x))},
    "2" = {overlap_obj <- function(x) length(unique(x)) == 1},
    "3" = {overlap_obj <- function(x) -length(unique(x))},
    "4" = {overlap_obj <- function(x) max(tabulate(x)) +
      2*(length(unique(x)) == 1)},
    "5" = {overlap_obj <- function(x) max(tabulate(x)) +
      -length(unique(x))}
  )

  all_comms <- list()
  all_comms[[1]] <- comms

  # store original overlap
  overlap <- vector("numeric", nrow(comms)+1)
  overlap[1] <- sum(apply(comms, 1, overlap_obj))

  # # swap the mode of each row for nonmatching labels, iterate over rows
  # # keep new ordering if objective is higher, otherwise discard
  # init <- sample(1:nrow(comms), match.iter, replace = TRUE)
  # for (i in 1:match.iter) {
  #   comms_new <- comms
  #   r <- init[i]
  #   x <- comms[r,]
  #   new <- as.numeric(names(which.max(table(x)))) # get mode
  #   repl <- which(x != new)
  #   for (j in repl) {
  #     reord <- 1:max(comms_new[,j]) # get reordering vector
  #     reord[x[j]] <- new # swap nonmatching with the most frequent
  #     reord[new] <- x[j]
  #     comms_new[,j] <- c(reord)[comms_new[,j]] # recode labels based on reord
  #   }
  #   overlap[[i+1]] <- sum(apply(comms_new, 1, overlap_obj))
  #   # replace comms if obj increases
  #   if (overlap[[i+1]] > overlap[[i]]) {comms <- comms_new}
  #   all_comms[[i+1]] <- comms_new # store even if not kept
  # }

  # replace on each iteration regardless
  init <- sample(1:nrow(comms))
  for (i in 1:nrow(comms)) {
    comms_new <- comms
    r <- init[i]
    x <- comms[r,]
    new <- as.numeric(names(which.max(table(x)))) # get mode
    repl <- which(x != new)
    for (j in repl) {
      reord <- 1:max(comms_new[,j]) # get reordering vector
      reord[x[j]] <- new # swap nonmatching with the most frequent
      reord[new] <- x[j]
      comms_new[,j] <- c(reord)[comms_new[,j]] # recode labels based on reord
    }
    overlap[[i+1]] <- sum(apply(comms_new, 1, overlap_obj))
    # replace no matter what
    comms <- comms_new
    all_comms[[i+1]] <- comms_new # store even if not kept
  }
  comms <- all_comms[[which.max(overlap)]]

  if(!is.null(manual.labs)) {comms[] <- sapply(comms,
                                               function(x) c(manual.labs)[x])}
  
  # order ROIs based first on most matching labels
  # then on number of times the mode in each row appears in whole data
  # then on each column
  agree <- apply(comms, 1, function(x) -max(table(x)))
  n_mode <- apply(comms, 1, function(x) {
    mode <- as.numeric(names(which.max(table(x))))
    -sum(c(comms) == mode)
  })
  comms_cols <- as.data.frame(comms)

  if (!is.null(atlas)) {
    comms <- cbind(atlas, comms)
    colnames(comms)[1] <- atlas.lab
    comms <- comms[do.call(order, c(list(atlas, agree, n_mode), comms_cols)),]
  } else {
    comms <- comms[do.call(order, c(list(agree, n_mode), comms_cols)),]
  }

  comms_m <- melt(comms)
  
  if (is.null(colors)) {
    ncolors <- length(unique(c(comms)))
    colors <- colorRampPalette(brewer.pal(8, "Set1"))(ncolors)
  }
  ggplot(data = comms_m) + geom_tile(aes(x = Var2, y = Var1, fill = as.factor(value))) +
    guides(fill=guide_legend(ncol=1)) +
    labs(x = "", y = "ROI", fill = "Community") +
    scale_fill_manual(values = colors,
      breaks = sort(unique(atlas))) +
    scale_y_reverse() +
    theme_classic() +
    theme(axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank())
}

getPerms <- function(x) {
  if (length(x) == 1) {
    return(x)
  }
  else {
    res <- matrix(nrow = 0, ncol = length(x))
    for (i in seq_along(x)) {
      res <- rbind(res, cbind(x[i], Recall(x[-i])))
    }
    return(res)
  }
}

# function for getting pairwise distances between site-specific clusterings
grp_comp <- function(matlab_res, grps, metric) {
  lapply(matlab_res, function(y) {
    x <- y[[1]]
    comp <- matrix(0, nrow(x), nrow(x),
                   dimnames = list(levels(grps), levels(grps)))
    for (i in 1:nrow(x)) {
      for (j in 1:nrow(x)) {
        if (i == j) {
          comp[i,j] <- 0
        } else {
          comp[i,j] = igraph::compare(x[i,], x[j,], metric)
        }
      }
    }
    comp
  })
}
