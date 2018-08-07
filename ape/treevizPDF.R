treevizPDF <- function(x, ..., height = NULL, width = 10, new = FALSE)
{
  if (!new) {
    if (!exists("last_treeviz", envir = .PlotPhyloEnv))
      new <- TRUE
    else
      tf <- get("last_treeviz", envir = .PlotPhyloEnv)$PDFfilename
  }
  if (new) {
    tf <- paste(tempfile(), "pdf", sep = ".")
    assign("last_treeviz", list(PDFfilename = tf),
           envir = .PlotPhyloEnv)
  }
  
  n <- Ntip(x)
  if (is.null(height)) height <- n * 0.19
  ## 0.19 is strheight("x", "inches") * 1.5
  pdf(tf, height = height, width = width)
  
  ## to avoid wide empty areas above and below the tree:
  par(yaxs = "i")
  plot.phylo(x, y.lim = c(-1, n + 1), ...)
  
  dev.off()
  system(paste(options()$pdfviewer, tf), wait = FALSE)
}

treevizGL <- function(x, edge.color = "orange", tip.color = "black")
{
  if (!require(rgl))
    stop("package 'rgl' not installed")
  
  rgl.segments <- function(x0, y0, x1, y1, z0 = NULL, z1 = NULL) {
    N <- length(x0)
    if (is.null(z0)) z0 <- rep(0, N)
    if (is.null(z1)) z1 <- rep(0, N)
    x <- as.vector(rbind(x0, x1))
    y <- as.vector(rbind(y0, y1))
    z <- as.vector(rbind(z0, z1))
    segments3d(x, y, z, col = edge.color)
  }
  
  rgl.cladogram.plot <- function(edge, xx, yy) {
    rgl.segments(xx[edge[, 1]], yy[edge[, 1]],
                 xx[edge[, 2]], yy[edge[, 2]])
    rgl.texts(xx[1:n], yy[1:n], rep(0, n), x$tip.label,
              adj = 0, col = tip.color)
  }
  
  n <- length(x$tip.label)
  m <- Nnode(x)
  N <- Nedge(x)
  
  yy <- numeric(n + m)
  TIPS <- x$edge[x$edge[, 2] <= n, 2]
  yy[TIPS] <- 1:n
  
  x <- reorder(x, order = "pruningwise")
  ans <- .C("node_height_clado", as.integer(n),
            as.integer(m), as.integer(x$edge[, 1]),
            as.integer(x$edge[, 2]), as.integer(N),
            double(n + m), as.double(yy),
            DUP = FALSE, PACKAGE = "ape")
  xx <- ans[[6]] - 1
  yy <- ans[[7]]
  
  rgl.open()
  rgl.cladogram.plot(x$edge, xx, yy)
}
