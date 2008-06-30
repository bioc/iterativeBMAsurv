\name{imageplot.iterate.bma.surv}
\alias{imageplot.iterate.bma.surv}
\title{An image plot visualization tool}
\description{Create a visualization of the models and variables selected by the
             iterative BMA algorithm.}
             
\usage{imageplot.iterate.bma.surv (bicreg.out, color="default", ...)}

\arguments{
\item{bicreg.out}{An object of type 'bicreg', 'bic.glm' or 'bic.surv'}
\item{color}{The color of the plot. The value 'default' uses the 
             current default R color scheme for image. The value 
             'blackandwhite' produces a black and white image.}
\item{...}{Other parameters to be passed to the image and axis functions.}
}

\details{This function is a modification of the \code{imageplot.bma}
         function from the \code{BMA} package. The difference is that
	 variables (genes) with \code{probne0} equal to 0 are removed
	 before plotting. The arguments of this function are identical
	 to those in \code{imageplot.bma}.}

\value{An heatmap-style image, with the BMA selected variables on the vertical
       axis, and the BMA selected models on the horizontal axis. The variables
       (genes) are sorted in descreasing order of the posterior probability 
       that the variable is not equal to 0 (\code{probne0}) from top to 
       bottom.  The models are sorted in descreasing order of the
       model posterior probability (\code{postprob}) from left to right.
       }

\references{
Annest, A., Yeung, K.Y., Bumgarner, R.E., and Raftery, A.E. (2008).
Iterative Bayesian Model Averaging for Survival Analysis.
Manuscript in Progress.

Clyde, M. (1999) 
Bayesian Model Averaging and Model Search Strategies (with discussion). In Bayesian Statistics 6. J.M. Bernardo, A.P. Dawid, J.O. Berger, and A.F.M. Smith eds. Oxford University Press, pages 157-185.

Yeung, K.Y., Bumgarner, R.E. and Raftery, A.E. (2005) 
Bayesian Model Averaging: Development of an improved multi-class, gene selection and classification tool for microarray data. 
Bioinformatics 21: 2394-2402.
}
\note{The \code{BMA} package is required.}

\seealso{\code{\link{iterateBMAsurv.train.wrapper}},
         \code{\link{iterateBMAsurv.train.predict.assess}},
         \code{\link{trainData}},
	 \code{\link{trainSurv}}, 
         \code{\link{trainCens}} 
}

\examples{
library (BMA)
library (iterativeBMAsurv)
data(trainData)
data(trainSurv)
data(trainCens)

## Training phase: select relevant genes
## Assumes the training data is in sorted order with the desired number of genes
ret.bic.surv <- iterateBMAsurv.train.wrapper (x=trainData, surv.time=trainSurv, cens.vec=trainCens)

## Produce an image plot to visualize the selected genes and models
imageplot.iterate.bma.surv (ret.bic.surv$obj)

}
\keyword{survival}



