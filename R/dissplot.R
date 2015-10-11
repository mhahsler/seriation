#######################################################################
# seriation - Infrastructure for seriation
# Copyrigth (C) 2011 Michael Hahsler, Christian Buchta and Kurt Hornik
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.



## Cluster visualization by proximity matrix shading


## interface
dissplot <- function(x, labels = NULL, method = "ARSA",
  control = NULL, options = NULL, ...) {
  
  ## add ... to options
  options <- c(options, list(...))

  
  
  
  ## make x dist
  if(!inherits(x, "dist")) {
    if(is.matrix(x) && isSymmetric(x)) x <- as.dist(x)
    else
      stop("Argument 'x' cannot safely be coerced to class 'dist'.")
  }
  
  res <- .arrange_dissimilarity_matrix(x, labels = labels,
    method = method, control = control)
  
  ## supress plot?
  plot <- if(is.null(options$plot)) TRUE else is.null(options$plot)
  if(plot) plot(res, options, gp = options$gp)
  
  invisible(res)
}


## work horse
.arrange_dissimilarity_matrix <- function(x, labels = NULL, method = NULL, 
  control = NULL) {
  
  ## x is already of class dist
  dim <- attr(x, "Size")
  diss_measure <- attr(x, "method")
  
  ## check labels
  if(!is.null(labels) && length(labels) != dim) 
    stop("Number of labels in 'labels' does not match dimensions of 'x'.")
  
  m <- method
  
  
  ## set everything to NULL first
  order               <- NULL
  k                   <- NULL             # number of clusters
  sil                 <- NULL
  avgSil              <- NULL
  labels_unique       <- NULL
  cluster_dissimilarities <- NULL
  ## method$a means method$ aggregation (default is avg)
  aggregation	<- "avg"
  if(class(method) == "list" && !is.null(method$a)) aggregation <- method$a
  
  if(class(method) != "list") method <- 
    list(inter_cluster = m, intra_cluster = m)
  m <- pmatch(names(method), c("inter_cluster", "intra_cluster", "aggregation"))
  if(any(is.na(m))) stop("Unknown method component. Use 'inter_cluster', 'intra_cluster' and 'aggregation'.")
  names(method) <- c("inter_cluster", "intra_cluster", "aggregation")[m]
  
  if(class(control[[1]]) != "list"){
    control <- list(inter_cluster = control, intra_cluster = control)
  }
  
  if(!is.null(method$inter_cluster) 
    && is.na(method$inter_cluster)) {
    ## no setiation
    if(!is.null(labels)) {
      ## do coarse seriation
      order <- order(labels)
      k <- length(unique(labels))
      ## calculate cluster_dissimilarities for later
      cluster_dissimilarities <- .cluster_dissimilarity(x, labels, 
        aggregation)
      aggregation <- attr(cluster_dissimilarities, "method")
      ## calculate silhouette values for later use
      sil <- cluster::silhouette(labels, x)
      
    }
    ## else keep the matrix as is -- do not reorder
    
  }else if(is.null(labels)) {
    ## reorder whole matrix if no labels are given
    order <- seriate(x, method = method$inter_cluster, 
      control = control$inter)[[1]] 
    
    method$inter_cluster <- if(!is.null(attr(order, "method"))) 
      attr(order, "method") else method$inter_cluster
    
    order <- get_order(order)
    
  }else{
    ## reorder clusters for given labels
    ## get number of clusters k
    k <- length(unique(labels))
    
    ## reorder with average pairwise dissimilarites between clusters
    cluster_dissimilarities <- .cluster_dissimilarity(x, labels, 
      aggregation)
    aggregation <- attr(cluster_dissimilarities, "method")
    
    if(k>2) {
      cluster_order <- seriate(as.dist(cluster_dissimilarities), 
        method = method$inter_cluster, control = control$inter)[[1]]
      
      method$inter_cluster <- if(!is.null(attr(cluster_order, "method"))) 
        attr(cluster_order, "method") else method$inter_cluster
      
      cluster_order <- get_order(cluster_order)
    }else{
      cluster_order <- 1:k
    }
    
    ## calculate silhouette values for later use
    sil <- cluster::silhouette(labels, x)
    
    ## determine order for matrix from cluster order
    order <- c()
    
    if(!is.null(method$intra_cluster) && 
        is.na(method$intra_cluster)) {
      ## no intra cluster ordering
      for(i in 1 : k) {
        order <- c(order, which(labels == cluster_order[i]))
      }
      ##method$intra_cluster <- NA
      
    }else{
      ## intra cluster order
      
      for(i in 1 : k) {
        take <- which(labels == cluster_order[i])
        
        ## only reorder for >1 elements
        if(length(take) > 1) {
          
          if(is.character(method$intra_cluster) &&
              match(tolower(method$intra_cluster), 
                c("sil", "silhouette", "silhouette width"), 
                nomatch = 0) > 0) {
            intra_order <-  order(sil[take, "sil_width"], 
              decreasing = TRUE)
            
            method$intra_cluster <- "silhouette width"
          }else{
            ## we use .rearrange_dist instead of permute
            ## since we take only a subset!
            block <- .rearrange_dist(x, take)
            
            intra_order <- seriate(block, method = method$intra_cluster, 
              control = control$intra)[[1]]
            
            method$intra_cluster <- 
              if(!is.null(attr(intra_order, "method")))
                attr(intra_order, "method") else method$intra_cluster
            
            intra_order <- get_order(intra_order)
          }
          
          order <- c(order, take[intra_order])
          
        }else{
          order <- c(order, take)
        }
        
      }
    }
    
    
    ## reorder cluster_dissimilarities for later
    cluster_dissimilarities  <- 
      cluster_dissimilarities[cluster_order, cluster_order]
    
  }
  
  ## reorder matrix
  if(!is.null(order)) {
    x_reordered <- permute(x, order)
    labels <- labels[order]
  }
  else x_reordered <- x
  
  ## prepare for return value
  cluster_description <- NULL
  
  if(!is.null(labels)) {
    
    labels_unique   <-  unique(labels)
    
    ## reorder silhouettes
    sil <- sil[order,]
    
    ## calculate avg silhouettes
    avgSil <- sapply(labels_unique, function(x) 
      mean(sil[sil[,"cluster"]==x, "sil_width"])) 
    
    ## generate description
    cluster_description = data.frame(
      position        = c(1 : k),
      label           = labels_unique, 
      size            = tabulate(labels)[labels_unique],
      ## FIXME: this is not the average anymore!
      aggregated_dissimilarity = 
        diag(cluster_dissimilarities)[labels_unique],
      avg_silhouette_width = avgSil)
  }
  
  ## clean order from names, etc.
  attributes(order) <- NULL
  
  result <- list(
    x_reordered     = x_reordered, 
    labels          = labels, 
    seriation_methods          = method, 
    aggregation_method          = aggregation, 
    k               = k, 
    cluster_dissimilarities =  cluster_dissimilarities,
    sil             = sil,
    order           = order, 
    cluster_order   = labels_unique,
    diss_measure    = diss_measure,
    description     =  cluster_description)
  
  class(result) <- "reordered_cluster_dissimilarity_matrix" 
  invisible(result)
}



## plot for reordered_cluster_dissimilarity_matrix
plot.reordered_cluster_dissimilarity_matrix <- function(x, options = NULL, ...) {
  
  m       <- as.matrix(x$x_reordered)
  k       <- x$k
  dim     <- attr(x$x_reordered, "Size")
  labels  <- x$labels
  labels_unique <- unique(labels)
  
  ## default plot options
  options <- .get_parameters(options, list(
    cluster_labels = TRUE, 
    lines       = TRUE, 
    averages	  = c(FALSE, TRUE), ## (upper.tri, lower.tri); NA means white
    flip		    = FALSE,
    silhouettes = FALSE,
    col         = NULL, 
    threshold   = NULL,
    zlim        = NULL,
    key         = TRUE,
    main        = NULL,
    axes        = "auto",
    gp          = gpar(),
    gp_lines    = gpar(),
    gp_labels   = gpar(),
    newpage     = TRUE, 
    pop         = TRUE
  )) 
  
  if(is.null(options$col)) 
    options$col <- rev(.sequential_pal())
  else 
    options$col <- rev(options$col)
  
  i <- pmatch(options$axes, c("auto", "x", "y", "both", "none"))
  if(is.na(i)) stop("Illegal vaule for axes. Use: 'auto', 'x', 'y', 'both' or 'none'!")
  options$axes <- c("auto", "x", "y", "both", "none")[i]
  
  ## clear page
  if(options$newpage) grid.newpage()
  
  ## do we have silhouettes?
  if(is.null(x$sil)) options$silhouettes <- FALSE
  
  ## create panels with avg. dissimilarity
  
  ## blank out if NA
  if(is.na(options$averages[1])) m[upper.tri(m)] <- NA
  if(is.na(options$averages[2])) m[lower.tri(m)] <- NA
  options$averages[is.na(options$averages)] <- FALSE
  
  if(!is.null(x$cluster_dissimilarities) 
    && !is.null(labels)
    && any(options$averages)){
    
    for(i in 1 : k) {
      for( j in 1 : k) {
        
        ## check empty clusters
        if(is.na(labels_unique[i])) next
        if(is.na(labels_unique[j])) next
        
        ## upper panels stay unchanged
        if(i < j && options$averages[1]) { 
          m[labels == labels_unique[i], labels == labels_unique[j]] <-
            x$cluster_dissimilarities[i, j] 
        }
        
        ## do lower panels
        if(i > j && options$averages[2]) { 
          m[labels == labels_unique[i], labels == labels_unique[j]] <-
            x$cluster_dissimilarities[i, j] 
        }
        
        ## do diagonal
        if(i == j) {
          block <- m[labels == labels_unique[i], 
            labels == labels_unique[j]]
          
          
          if(options$averages[1]) { 
            block[upper.tri(block, diag = TRUE)] <- 
              x$cluster_dissimilarities[i, j]
            
            m[labels == labels_unique[i],
              labels == labels_unique[j]] <- block
          }
          
          if(options$averages[2]) { 
            block[lower.tri(block, diag = TRUE)] <- 
              x$cluster_dissimilarities[i, j]
            
            m[labels == labels_unique[i],
              labels == labels_unique[j]] <- block
          }
          
        }
      }
    }
  }
  
  if(options$flip){
    m <- m[,ncol(m):1]
  }
  
  if(!options$silhouettes) {
    pushViewport(viewport(layout = grid.layout(6, 3,
      widths = unit.c(
        unit(2, "lines"),                       # space
        unit(1, "snpc") - unit(7, "lines"),     # image
        unit(2, "lines")                        # space
      ),
      heights = unit.c(
        unit(2, "lines"),                       # title
        unit(1, "lines"),                       # space
        unit(1, "snpc") - unit(7, "lines"),     # image
        unit(1, "lines"),                       # space
        unit(1, "lines"),                       # colorkey
        unit(2, "lines")                        # space
      )
    ), gp=options$gp))
    
    main_vp     <- viewport(layout.pos.col = 2, layout.pos.row = 1, 
      name="main")
    image_vp    <- viewport(layout.pos.col = 2, layout.pos.row = 3)
    colorkey_vp <- viewport(layout.pos.col = 2, layout.pos.row = 5,
      name = "colorkey")
    
  }else{ ## with silhouettes
    pushViewport(viewport(layout = grid.layout(6, 5,
      widths = unit.c(
        unit(2, "lines"),                       # space
        unit(0.7, "snpc") - unit(2.5, "lines"), # image
        unit(1, "lines"),                       # space
        unit(0.3, "snpc") - unit(2.5, "lines"), # sil
        unit(2, "lines")                        # space
      ),
      heights = unit.c(
        unit(2, "lines"),                       # title
        unit(2, "lines"),                       # space
        unit(0.7, "snpc") - unit(2.5, "lines"), # image
        unit(1, "lines"),                       # space
        unit(1, "lines"),                       # colorkey
        unit(2, "lines")                        # space
      )
    ), gp = options$gp))
    
    main_vp     <- viewport(layout.pos.col = 2:4, layout.pos.row = 1,
      name="main")
    image_vp    <- viewport(layout.pos.col = 2,   layout.pos.row = 3)
    sil_vp  <- viewport(layout.pos.col = 4,   layout.pos.row = 3,
      name="sil")
    colorkey_vp <- viewport(layout.pos.col = 2,   layout.pos.row = 5,
      name="colorkey")
    
  }
  
  ## main
  pushViewport(main_vp)
  grid.text(options$main, gp = gpar(cex = 1.3, fontface = "bold"))
  upViewport(1)
  
  ## silhouette
  if(options$silhouettes) {
    ## get and reorder silhouettes
    s <- x$sil[,"sil_width"]
    
    pushViewport(sil_vp)
    .grid_barplot_horiz(s, xlab = "Silhouette width", 
      gp_bars = gpar(fill = "lightgrey", col = 0))
    upViewport(1)
    
  }
  
  ## image
  if(is.null(options$zlim)) options$zlim <- range(m, na.rm = TRUE)
  if(!is.null(options$threshold)) m[m > options$threshold] <- NA

  pushViewport(image_vp)
  .grid_image(m, col = options$col, zlim = options$zlim)
  
  ## add labels?
  if(options$axes == "auto" && nrow(m)>25) options$axes <- "none"
  if(options$axes != "none") {
    downViewport("image")
    #grid.text(colnames(m), y = unit(-1, "lines"),
    #  x=unit(1:ncol(m), "native"), rot=90, just="right")
    
    grid.text(rownames(m), x = unit(1, "npc") + unit(1, "lines"),
      y=unit(1:nrow(m), "native"), just="left", gp = options$gp_labels)
    upViewport(1)
  }
  
  upViewport(1)
  
  ## color key?  
  if(options$key) {
    pushViewport(colorkey_vp)
    .grid_colorkey(options$zlim, col = options$col, 
      threshold = options$threshold)
    upViewport(1)
  }
  
  ## plot cluster borders if we have labels and order
  if(!is.null(labels)) {
    
    labels_unique_y		<- labels_unique
    cluster_width_y		<- (tabulate(labels)[labels_unique])
    #cluster_cuts_y		<- cumsum(cluster_width_y) + 0.5
    cluster_cuts_y		<- cumsum(cluster_width_y)
    cluster_center_y	<- cluster_cuts_y - cluster_width_y / 2
    
    if(options$flip) {
      labels_unique_x	<- rev(labels_unique)
      cluster_width_x <- (tabulate(labels)[labels_unique_x])
      #	cluster_cuts_x  <- cumsum(cluster_width_x) + 0.5
      cluster_cuts_x  <- cumsum(cluster_width_x)
      cluster_center_x<- cluster_cuts_x - cluster_width_x / 2
    }else{
      labels_unique_x <- labels_unique_y
      cluster_width_x <- cluster_width_y
      cluster_cuts_x <- cluster_cuts_y
      cluster_center_x <- cluster_center_y
    }
    
    if(options$cluster_labels) {
      
      seekViewport("image")
      
      ## above the plot
      grid.text(labels_unique_x,
        x = cluster_center_x, 
        y = unit(1, "npc") + unit(1, "lines"),
        default.units = "native", gp = options$gp_labels)
      ## left of the plot
      grid.text(labels_unique_y, 
        x = unit(-1, "lines"),
        y = cluster_center_y,
        default.units = "native", gp = options$gp_labels)
      upViewport(2)
    }
    
    if(options$lines) {
      ## draw lines separating the clusters
      #cluster_cuts <- cluster_cuts[-length(cluster_cuts)]
      ## remove last line
      
      seekViewport("image")
      for(i in 1:(k-1)) {
        
        grid.lines(
          #x = c(0, dim), 
          x = c(0.5, dim + 0.5), 
          y = cluster_cuts_y[i]+.5, 
          default.units = "native", gp = options$gp_lines)
        
        grid.lines(
          x = cluster_cuts_x[i]+.5, 
          #y = c(0, dim), 
          y = c(0.5, dim + 0.5), 
          default.units = "native", gp = options$gp_lines)
        
      }
      upViewport(2)
      
    }
  }
  
  if (options$pop) popViewport(1)
  else upViewport(1)
}


## print for reordered_cluster_dissimilarity_matrix
print.reordered_cluster_dissimilarity_matrix <-
  function(x, ...)
  {
    d <- attr(x$x_reordered, "Size")
    k <- if(!is.null(x$k)) x$k else NA
    
    cat(gettextf("object of class '%s'\n", class(x)))
    cat("matrix dimensions:", d, "x", d, "\n")
    cat(gettextf("dissimilarity measure: '%s'\n", x$diss_measure))
    cat("number of clusters k:", k, "\n")
    if(!is.null(x$k)) {
      cat("\ncluster description\n")
      print(x$description)
    }
    
    cat("\n")
    cat("used seriation methods\n")
    cat(gettextf("inter-cluster: '%s'\n", x$seriation_methods$inter))
    cat(gettextf("intra-cluster: '%s'\n", x$seriation_methods$intra))
    
    cat("\n")
    cat(gettextf("dissimilarity aggregation method: '%s'\n", 
      x$aggregation_method))
    
    invisible(x)
  }

## inter and intra cluster dissimilarity matrix from 
## a dissimilarity matrix plus labels
.cluster_dissimilarity <- function(x, labels, method=c("avg", "min", "max",
  "Hausdorff")) {
  
  method <- match.arg(method)
  ## FIXME: Implement Hausdorff
  
  linkage <- if(method=="avg") mean
  else if(method=="min") min
  else if(method=="max") max
  else if(method=="Hausdorff") .hausdorff
  else stop("Unknown method.")
  
  if(class(x) != "matrix") x <- as.matrix(x)
  
  
  ## kill self-dissimilarities (which are always 0)
  diag(x) <- NA
  
  k           <- length(unique(labels))
  diss_matrix <- matrix(nrow = k, ncol = k)
  
  ## calculate avg. dissimilarity between clusters
  for(i in 1:k) {
    slice <- x[labels == i, , drop = FALSE]
    for(j in 1:i) {
      block <- slice[,labels == j, drop = FALSE]
      
      val <- linkage(block, na.rm = TRUE)
      
      ## fix for clusters of size 1
      if(is.nan(val)) val <- 0
      
      diss_matrix[i, j] <- val
      diss_matrix[j, i] <- val
    }
  }
  
  attr(diss_matrix, "method") <- method
  diss_matrix
}

## implement Hausdorff distance between two sets from a dissimilarity matrix
##d_H = max{sup_x\inX inf_y\inY d(x,y), sup_y\inY inf_x\inX d(x,y)} 
.hausdorff <- function(block, na.rm=TRUE) max(
  apply(block, MARGIN = 1, min, na.rm = na.rm),
  apply(block, MARGIN = 2, min, na.rm = na.rm))

