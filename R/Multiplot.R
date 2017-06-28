#' Multiple plot function
#'
#' @param ... ggplot objects
#' @param plotlist ggplot objects
#' @param file file
#' @param cols cols
#' @param layout layout (see description)
#'
#' @description ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
#'
#' - cols:   Number of columns in layout
#'
#' - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#'
#' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#'
#' @export
#'
#' @examples library(ggplot2)
#'
#' This example uses the ChickWeight dataset, which comes with ggplot2
#' First plot
#' p1 <- ggplot(ChickWeight, aes(x=Time, y=weight, colour=Diet, group=Chick)) +
#'  geom_line() +
#'  ggtitle("Growth curve for individual chicks")
#'
#' Second plot
#' p2 <- ggplot(ChickWeight, aes(x=Time, y=weight, colour=Diet)) +
#'  geom_point(alpha=.3) +
#'  geom_smooth(alpha=.2, size=1) +
#'  ggtitle("Fitted growth curve per diet")
#'
# Third plot
#' p3 <- ggplot(subset(ChickWeight, Time==21), aes(x=weight, colour=Diet)) +
#'  geom_density() +
#'  ggtitle("Final weight, by diet")
#'
#' Fourth plot
#' p4 <- ggplot(subset(ChickWeight, Time==21), aes(x=weight, fill=Diet)) +
#'  geom_histogram(colour="black", binwidth=50) +
#'  facet_grid(Diet ~ .) +
#'  ggtitle("Final weight, by diet") +
#'  theme(legend.position="none")        # No legend (redundant in this graph)
#'
#'  multiplot(p1, p2, p3, p4, cols=2)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' ggplot_missing
#'
#' @param x dataframe
#'
#' @return ranking variables with missing values
#'
#' @export
#'
ggplot_missing <- function(x){

  colSums(is.na(x)) %>%
    t %>% as.data.frame() %>%
    gather() %>%
    ggplot(aes(x = reorder(key,desc(value)), y = value / nrow(x) * 100)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
    labs(x = "Variables in Dataset",
         y = "observations / Rows %") +
    ggtitle("PERCENTAGE OF MISSING VALUES FOR EACH VARIABLE") +
    scale_y_continuous(limits = c(0,101), expand = c(0, 0))
}


# Quantile-comparison plots (J. Fox)

#' Quantile-comparison plots (J. Fox)
#'
#' @export
#'
qqPlot = function(object, ...) UseMethod("qqPlot", object)

qqPlot.default = function(x, distribution="norm", ..., ylab=deparse(substitute(x)),
                          xlab=paste(distribution, "quantiles"), main=NULL, las=par("las"),
                          envelope=.95,
                          col=palette()[1], col.lines=palette()[2], lwd=2, pch=1, cex=par("cex"),
                          line=c("quartiles", "robust", "none"),
                          labels = if(!is.null(names(x))) names(x) else seq(along=x),
                          id.method = "y",
                          id.n = if(id.method[1]=="identify") Inf else 0,
                          id.cex=1, id.col=palette()[1], id.location="lr", grid=TRUE)
{
  line <- match.arg(line)
  good <- !is.na(x)
  ord <- order(x[good])
  if (length(col) == length(x)) col <- col[good][ord]
  if (length(pch) == length(x)) pch <- pch[good][ord]
  if (length(cex) == length(x)) cex <- cex[good][ord]
  ord.x <- x[good][ord]
  ord.lab <- labels[good][ord]
  q.function <- eval(parse(text=paste("q", distribution, sep="")))
  d.function <- eval(parse(text=paste("d", distribution, sep="")))
  n <- length(ord.x)
  P <- ppoints(n)
  z <- q.function(P, ...)
  plot(z, ord.x, type="n", xlab=xlab, ylab=ylab, main=main, las=las)
  if(grid){
    grid(lty=1, equilogs=FALSE)
    box()}
  points(z, ord.x, col=col, pch=pch, cex=cex)
  if (line == "quartiles" || line == "none"){
    Q.x <- quantile(ord.x, c(.25,.75))
    Q.z <- q.function(c(.25,.75), ...)
    b <- (Q.x[2] - Q.x[1])/(Q.z[2] - Q.z[1])
    a <- Q.x[1] - b*Q.z[1]
    if (line == "quartiles") abline(a, b, col=col.lines, lwd=lwd)
  }
  if (line=="robust") {
    coef <- coef(rlm(ord.x ~ z))
    a <- coef[1]
    b <- coef[2]
    abline(a, b, col=col.lines, lwd=lwd)
  }
  conf <- if (envelope == FALSE) .95 else envelope
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (b/d.function(z, ...))*sqrt(P*(1 - P)/n)
  fit.value <- a + b*z
  upper <- fit.value + zz*SE
  lower <- fit.value - zz*SE
  if (envelope != FALSE) {
    lines(z, upper, lty=2, lwd=lwd, col=col.lines)
    lines(z, lower, lty=2, lwd=lwd, col=col.lines)
  }
  showLabels(z, ord.x, labels=ord.lab,
             id.method = id.method, id.n = id.n, id.cex=id.cex, id.col=id.col, id.location=id.location)
}

qqPlot.lm <- function(x, xlab=paste(distribution, "Quantiles"),
                      ylab=paste("Studentized Residuals(", deparse(substitute(x)), ")", sep=""), main=NULL,
                      distribution=c("t", "norm"), line=c("robust", "quartiles", "none"), las=par("las"),
                      simulate=TRUE, envelope=.95,  reps=100,
                      col=palette()[1], col.lines=palette()[2], lwd=2, pch=1, cex=par("cex"),
                      labels, id.method = "y",
                      id.n = if(id.method[1]=="identify") Inf else 0, id.cex=1,
                      id.col=palette()[1], id.location="lr", grid=TRUE, ...){
  result <- NULL
  distribution <- match.arg(distribution)
  line <- match.arg(line)
  rstudent <- rstudent(x)
  if (missing(labels)) labels <- names(rstudent)
  good <- !is.na(rstudent)
  rstudent <- rstudent[good]
  labels <- labels[good]
  sumry <- summary.lm(x)
  res.df <- sumry$df[2]
  if(!simulate)
    result <- qqPlot(rstudent, distribution=if (distribution == "t") "t" else "norm", df=res.df-1, line=line,
                     main=main, xlab=xlab, ylab=ylab, las=las, envelope=envelope,
                     col=col, col.lines=col.lines, lwd=lwd, pch=pch, cex=cex,
                     labels=labels, id.method=id.method, id.n=id.n, id.cex=id.cex,
                     id.col=id.col, id.location="lr", ...)
  else {
    n <- length(rstudent)
    ord <- order(rstudent)
    ord.x <- rstudent[ord]
    ord.lab <- labels[ord]
    P <- ppoints(n)
    z <- if (distribution == 't') qt(P, df=res.df-1) else qnorm(P)
    plot(z, ord.x, type="n", xlab=xlab, ylab=ylab, main=main, las=las)
    if(grid) grid(lty=1, equilogs=FALSE)
    points(z, ord.x, pch=pch, col=col, cex=cex)
    yhat <- na.omit(fitted.values(x))
    S <- sumry$sigma
    Y <- matrix(yhat, n, reps) + matrix(rnorm(n*reps, sd=S), n, reps)
    X <- model.matrix(x)
    rstud <- apply(rstudent(lm(Y ~ X - 1)), 2, sort)
    lower <- apply(rstud, 1, quantile, prob=(1 - envelope)/2)
    upper <- apply(rstud, 1, quantile, prob=(1 + envelope)/2)
    lines(z, upper, lty=2, lwd=lwd, col=col.lines)
    lines(z, lower, lty=2, lwd=lwd, col=col.lines)
    if (line == "quartiles"){
      Q.x <- quantile(rstudent, c(.25,.75))
      Q.z <- if (distribution == 't') qt(c(.25,.75), df=res.df - 1) else qnorm(c(.25,.75))
      b <- (Q.x[2] - Q.x[1])/(Q.z[2] - Q.z[1])
      a <- Q.x[1] - b*Q.z[1]
      abline(a, b, col=col.lines, lwd=lwd)
    }
    if (line=="robust"){
      coef <- coefficients(rlm(ord.x~z))
      a <- coef[1]
      b <- coef[2]
      abline(a, b, col=col.lines, lwd=lwd)
    }
    result <- showLabels(z, ord.x,labels=ord.lab,
                         id.method = id.method, id.n = id.n, id.cex=id.cex, id.col=id.col, id.location=id.location)
  }
  if (length(result) == 0) invisible(result) else if (is.numeric(result)) sort(result) else result
}

qqPlot.glm <- function(x, ...) stop("QQ plot for studentized residuals not available for glm")

showLabels <- function(x, y, labels=NULL, id.method="identify",
                       id.n = length(x), id.cex=1, id.col=palette()[1], id.location="lr", ...) {
  res <- NULL
  id.method <- if(is.list(id.method)) id.method else list(id.method)
  for (meth in id.method)
    res <- c(res, showLabels1(x, y, labels, meth, id.n, id.cex,
                              id.col, id.location, ...))
  return(if(is.null(res)) invisible(res) else res)
}

showLabels1 <- function(x, y, labels=NULL, id.method="identify",
                        id.n = length(x), id.cex=1, id.col=palette()[1], id.location="lr", all=NULL, ...) {
  # If labels are NULL, try to get the labels from x:
  if (is.null(labels)) labels <- names(x)
  if (is.null(labels)) labels <- paste(seq_along(x))
  if (is.null(id.col)) id.col <- palette()[1]
  if (is.null(id.location)) id.location <- "lr"
  id.location <- match.arg(id.location, c("lr", "ab"))
  # logged-axes?
  log.x <- par("xlog")
  log.y <- par("ylog")
  # id.method can be any of the following:
  #    --- a list of row numbers
  #    --- a list of labels
  #    --- a vector of n numbers
  #    --- a text string:  'identify', 'x', 'y', 'mahal', 'r'
  idmeth <- pmatch(id.method[1], c("x", "y", "mahal", "identify", "r"))
  if(!is.na(idmeth)) idmeth <- c("x", "y", "mahal", "identify", "r")[idmeth]
  # if idmeth is NA, then it must be <= n numbers or row names
  id.var <- NULL
  if(is.na(idmeth)){
    if(is.null(all)) all <- list(labels=labels, subs=rep(TRUE, length(labels)))
    names(all$labels) <- all$labels
    if(length(id.method) >= length(x)){
      id.var <- id.method[which(all$subs)]
      id.n <- min(id.n, length(id.var))
    }
    else {
      id.var <- rep(0, length(x))
      names(id.var) <- labels
      inSubset <- all$labels[all$subs] %in% all$labels[id.method]
      id.var[inSubset] <- 1
      id.n <- sum(inSubset)
    }
  }
  else {
    # use identify?
    if(idmeth == "identify"){
      result <- labels[identify(x, y, labels, n=length(x), cex=id.cex,
                                col=id.col, ...)]
      if(length(result) > 0) return(unique(result)) else return(NULL)
    }
    # missing values need to be removed
    ismissing <- is.na(x) | is.na(y) | is.na(labels)
    if( any(ismissing) ) {
      x <- x[!ismissing]
      y <- y[!ismissing]
      labels <- labels[!ismissing]
    }
    # other methods:
    id.var <- switch(id.method,
                     x = if(log.x==TRUE)
                       suppressWarnings(if(all(x) > 0)
                         abs(log(x) - mean(log(x))) else
                           return(invisible(NULL)))  else
                             abs(x - mean(x)),
                     y = if(log.y==TRUE)
                       suppressWarnings(if(all(y) > 0)
                         abs(log(y) - mean(log(y))) else
                           return(invisible(NULL)))  else
                             abs(y - mean(y)),
                     r = if(log.y==TRUE)
                       suppressWarnings(if(all(y) > 0)
                         abs(log(y)) else
                           return(invisible(NULL)))  else
                             abs(y),
                     mahal = if(log.x == TRUE & log.y == TRUE) {
                       suppressWarnings(if(all(x) > 0 & all(y) > 0)
                         rowSums( qr.Q(qr(cbind(1, log(x), log(y))))^2 ) else
                           return(invisible(NULL))) } else {
                             if(log.x == TRUE) {
                               suppressWarnings(if(all(x) > 0 )
                                 rowSums( qr.Q(qr(cbind(1, log(x), y)))^2 ) else
                                   return(invisible(NULL))) } else {
                                     if(log.y == TRUE) {
                                       suppressWarnings(if(all(y) > 0 )
                                         rowSums( qr.Q(qr(cbind(1, x, log(y))))^2 ) else
                                           return(invisible(NULL)))  } else {
                                             rowSums( qr.Q(qr(cbind(1, x, y)))^2 ) }}})
  }
  # require id.n positive
  if(id.n <= 0L) return(invisible(NULL))
  # criterion
  ind <-  order(id.var, decreasing=TRUE)[1L:min(length(id.var), id.n)]
  # position, now depends on id.location (as of 5/16/2016)
  if(id.location %in% c("lr", "l", "r")){
    mid <- mean(if(par("xlog")==TRUE) 10^(par("usr")[1:2]) else
      par("usr")[1:2])
    labpos <- c(4,2)[1+as.numeric(x > mid)]
  } else {
    mid <- mean(if(par("ylog")==TRUE) 10^(par("usr")[3:4]) else
      par("usr")[3:4])
    labpos <- c(3,1)[1+as.numeric(y > mid)]
  }
  # print
  for (i in ind) {
    text(x[i], y[i], labels[i], cex = id.cex, xpd = TRUE,
         col = id.col, pos = labpos[i], offset = 0.25, ...)}
  names(ind) <- labels[ind]
  result <- ind
  if (length(result) == 0) return(NULL) else return(result)
}


qqPlot.saddlepoint = function(object, ...) {

  warning("Quantile plot in saddlepoint class is experimental")

  x = object@data[order(object@data)]
  n = length(x)
  p = seq(0, 1, length.out = n)

  z = unname(quantile(object@transform, p))

  plot(z, x, type = "n", xlab = "tranformation", ylab = "data", main = NULL, las = par("las"))
  grid(lty = 1, equilogs = FALSE)
  box()
  points(z, x, col= palette()[1], pch= 1, cex= par("cex"))

    Q.x = quantile(x, c(.25,.75))
    Q.z =  quantile(z, c(.25,.75))
    b = (Q.x[2] - Q.x[1])/(Q.z[2] - Q.z[1])
    a = Q.x[1] - b*Q.z[1]
    abline(a, b, col = palette()[2], lwd = 2)


  conf = 0.95
  zz <- qnorm(1 - (1 - conf)/2)

  # Funcion de densidad de z
  SE <- (b / dnorm(z, mean(z), sd = sd(z))) * sqrt(p * (1 - p)/n)

  fit.value <- a + b*z
  upper <- fit.value + zz * SE
  lower <- fit.value - zz * SE
    lines(z, upper, lty = 2, lwd = 2, col = palette()[2])
    lines(z, lower, lty = 2, lwd = 2, col = palette()[2])
  showLabels(z, x, labels= palette()[2],
             id.method = "y", id.n = 0, id.cex = 1, id.col = palette()[1], id.location = "lr")
}

