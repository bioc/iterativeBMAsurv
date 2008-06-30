
## Goal:      To plot image plot without variables with probne0 == 0
##            and variables sorted by probne0 in descending order
## Arguments: Same as imageplot.bma
## Calls:     imageplot.bma

imageplot.iterate.bma.surv <- function (bicreg.out, color="default", ...) {
   
   ## Get genes with non-zero probne0
   selected.genes <- bicreg.out$namesx [bicreg.out$probne0 > 0]
   selected.probne0 <- bicreg.out$probne0 [bicreg.out$probne0 > 0]
   
   ## Sort by probne0
   sorted.vec <- sort (selected.probne0, decreasing = TRUE, index = TRUE)
   sorted.genes <- selected.genes [sorted.vec$ix]

   ## Modify components in bicreg.out: namesx, probne0, n.vars, which
   bicreg.mod.out <- bicreg.out
   bicreg.mod.out$namesx <- sorted.genes
   bicreg.mod.out$probne0 <- sorted.vec$x
   bicreg.mod.out$n.vars <- length (sorted.genes)
   
   ## DO NOT sort the columns in which
   bicreg.mod.out$which <- bicreg.out$which [, bicreg.out$probne0 > 0]
   bicreg.mod.out$which <- bicreg.mod.out$which[, sorted.vec$ix]

   ## Calls imageplot.bma with the modified argument
   imageplot.bma.mod (bicreg.mod.out, color=color)
}


## Code copied from imageplot.bma from the BMA package
## with a few modifications: margin of plot and title of plot

imageplot.bma.mod <- function (bicreg.out, color = "default", ...)  {
    
    keep.mar <- par()$mar
    par(mar = c(5, 8, 4, 2) + 0.1)
    which.out <- bicreg.out$which
    nvar <- ncol(which.out)
    nmodel <- nrow(which.out)
    par(las = 1)
    if (color == "default") {
        image(c(0, cumsum(bicreg.out$postprob)), 1:nvar, -which.out[1:nmodel, 
            nvar:1, drop = FALSE], xlab = "Model #", ylab = "", 
            xaxt = "n", yaxt = "n", xlim = c(0, 1), main = "Models selected by iterativeBMAsurv", 
            ...)
    }
    if (color == "blackandwhite") {
        image(cumsum(bicreg.out$postprob), 1:nvar, -which.out[1:nmodel, 
            nvar:1, drop = FALSE], xlab = "Model #", ylab = "", 
            xaxt = "n", yaxt = "n", col = c("black", "white"), 
            main = "Models selected by BMA", ...)
    }
    xat <- (cumsum(bicreg.out$postprob) + c(0, cumsum(bicreg.out$postprob[-nmodel])))/2
    axis(1, at = xat, labels = 1:nmodel, ...)
    axis(2, at = 1:nvar, labels = rev(bicreg.out$namesx), ...)
    par(mar = keep.mar)
}

