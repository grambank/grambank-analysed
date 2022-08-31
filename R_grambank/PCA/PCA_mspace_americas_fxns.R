
##################################
######## SakamotoVenditti2018  ###
##################################


nodeAges <- function(tree, tip.age=NULL, min.age=NULL){
  pairdist <- dist.nodes(tree)[,length(tree$tip.label)+1]
  nodedist <- pairdist[-c(1:length(tree$tip.label))]
  termdist <- pairdist[-c((length(tree$tip.label)+1):length(pairdist))]
  distMin <- min(termdist)
  if(is.null(tip.age)){ages <- max(pairdist) - pairdist
  if(is.null(min.age)){ages <- ages}else{ages <- ages + min.age}
  tip.age <- ages[1:length(tree$tip.label)]
  nodeage <- ages[-c(1:length(tree$tip.label))]
  }else{nodeage <- (max(tip.age) + distMin) - nodedist}
  allages <- list(tip.age=tip.age, nodeage=nodeage)
  return(allages)
}


# colours
vector_col<-function(branchestocol, logcolor=F){
  dpal <- get("palettes", envir = BAMMtools:::.colorEnv)
  palfun <- colorRampPalette(dpal$RdYlBu, space="Lab")
  focal.breaks <- assignColorBreaks(branchestocol, NCOLORS=64, logcolor=logcolor, method="linear")
  
  focal.pal <- palfun(length(focal.breaks))
  focal.pal <- magma(length(focal.breaks))
  # focal.pal <- rev(viridis(length(focal.breaks)))
  # focal.pal <- ocean.dense(length(focal.breaks))
  
  cuts <- as.numeric(cut(branchestocol, breaks = focal.breaks, include.lowest=T))
  cols <- focal.pal[cuts]
  cols[is.na(cols)] <- focal.pal[length(focal.pal)]
  
  toret<-list(cols, focal.pal)
  names(toret)<-c("cols", "focalpal")
  return(toret)
}

# mod imageplot
myimageplot<-function (..., add = FALSE, breaks = NULL, nlevel = 64, col = NULL, 
                       horizontal = FALSE, legend.shrink = 0.9, legend.width = 1.2, 
                       legend.mar = ifelse(horizontal, 3.1, 5.1), legend.lab = NULL, 
                       legend.line = 2, graphics.reset = FALSE, bigplot = NULL, 
                       smallplot = NULL, legend.only = FALSE, lab.breaks = NULL, 
                       axis.args = NULL, legend.args = NULL, legend.cex = 1, midpoint = FALSE, 
                       border = NA, lwd = 1, verbose = FALSE) 
{
  old.par <- par(no.readonly = TRUE)
  if (is.null(col)) {
    col <- tim.colors(nlevel)
  }
  else {
    nlevel <- length(col)
  }
  info <- imagePlotInfo(..., breaks = breaks, nlevel = nlevel)
  breaks <- info$breaks
  if (verbose) {
    print(info)
  }
  if (add) {
    big.plot <- old.par$plt
  }
  if (legend.only) {
    graphics.reset <- TRUE
  }
  if (is.null(legend.mar)) {
    legend.mar <- ifelse(horizontal, 3.1, 5.1)
  }
  temp <- imageplot.setup(add = add, legend.shrink = legend.shrink, 
                          legend.width = legend.width, legend.mar = legend.mar, 
                          horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
  smallplot <- temp$smallplot
  bigplot <- temp$bigplot
  if (!legend.only) {
    if (!add) {
      par(plt = bigplot)
    }
    if (!info$poly.grid) {
      image(..., breaks = breaks, add = add, col = col)
    }
    else {
      poly.image(..., add = add, breaks = breaks, col = col, 
                 midpoint = midpoint, border = border, lwd.poly = lwd)
    }
    big.par <- par(no.readonly = TRUE)
  }
  if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
    par(old.par)
    stop("plot region too small to add legend\n")
  }
  ix <- 1:2
  iy <- breaks
  nBreaks <- length(breaks)
  midpoints <- (breaks[1:(nBreaks - 1)] + breaks[2:nBreaks])/2
  iz <- matrix(midpoints, nrow = 1, ncol = length(midpoints))
  if (verbose) {
    print(breaks)
    print(midpoints)
    print(ix)
    print(iy)
    print(iz)
    print(col)
  }
  par(new = TRUE, pty = "m", plt = smallplot, err = -1)
  if (!horizontal) {
    image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
          ylab = "", col = col, breaks = breaks,col.axis="white",col.lab="white")
  }
  else {
    image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
          ylab = "", col = col, breaks = breaks, col.axis="white",col.lab="white")
  }
  if (!is.null(lab.breaks)) {
    axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
                        mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2), 
                        at = breaks, labels = lab.breaks), axis.args,col.axis="white",col.lab="white")
  }
  else {
    axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
                        mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)), 
                   axis.args,col.axis="white",col.lab="white")
  }
  do.call("axis", axis.args)
  if (!is.null(legend.lab)) {
    legend.args <- list(text = legend.lab, side = ifelse(horizontal, 
                                                         1, 4), line = legend.line, cex = legend.cex, col="white")
  }
  if (!is.null(legend.args)) {
    do.call(mtext, legend.args)
  }
  mfg.save <- par()$mfg
  if (graphics.reset | add) {
    par(old.par)
    par(mfg = mfg.save, new = FALSE)
    invisible()
  }
  else {
    par(big.par)
    par(plt = big.par$plt, xpd = FALSE)
    par(mfg = mfg.save, new = FALSE)
    invisible()
  }
}


