\name{predictBicSurv}
\alias{predictBicSurv}
\title{Predicted patient risk scores from iterative Bayesian Model Averaging}
\description{This function predicts the risk scores for patient
             samples in the test set.}
             
\usage{predictBicSurv(newdata.vec, postprob.vec, mle.mat)}

\arguments{
\item{newdata.vec}{A vector consisting of the data from a test sample.}
\item{postprob.vec}{A vector consisting of the posterior probability
                    of each BMA selected model.}
\item{mle.mat}{A matrix with one row per model and one column per variable giving 
               the maximum likelihood estimate of each coefficient for each 
               \code{bic.surv} selected model.}
}

\details{The function begins by computing the risk score of each model k in 
         the selected set of models M. The risk score for a model k = sum
         (coefficient in model k * corresponding expression level in newdata.vec). 
         The function then predicts a patient risk score by summing the product 
         of the posterior probability of model k and the risk score of model k 
         over all models in M. In other words, predicted patient risk score = 
         sum (postprob model k * risk score model k). This function is called
         by \code{iterateBMAsurv.train.predict.assess}.}

\value{A real number representing the predicted risk score of a given 
       patient sample.}

\references{
Annest, A., Yeung, K.Y., Bumgarner, R.E., and Raftery, A.E. (2008).
Iterative Bayesian Model Averaging for Survival Analysis.
Manuscript in Progress.

Raftery, A.E. (1995). 
Bayesian model selection in social research (with Discussion). Sociological Methodology 1995 (Peter V. Marsden, ed.), pp. 111-196, Cambridge, Mass.: Blackwells.

Volinsky, C., Madigan, D., Raftery, A., and Kronmal, R. (1997)
Bayesian Model Averaging in Proprtional Hazard Models: Assessing the Risk of a Stroke. 
Applied Statistics 46: 433-448.

Yeung, K.Y., Bumgarner, R.E. and Raftery, A.E. (2005) 
Bayesian Model Averaging: Development of an improved multi-class, gene selection and classification tool for microarray data. 
Bioinformatics 21: 2394-2402.
}
\note{}

\seealso{\code{\link{iterateBMAsurv.train.predict.assess}},
         \code{\link{iterateBMAsurv.train.wrapper}},
         \code{\link{predictiveAssessCategory}},
         \code{\link{trainData}},
         \code{\link{trainSurv}},
         \code{\link{trainCens}},
         \code{\link{testData}}
}

\examples{
library (BMA)
library (iterativeBMAsurv)
data(trainData)
data(trainSurv)
data(trainCens)
data(testData)

## Training phase: select relevant genes. Assume the training data is sorted 
## and includes the desired number of top-ranked genes.
ret.list <- iterateBMAsurv.train.wrapper (x=trainData, surv.time=trainSurv, cens.vec=trainCens, nbest=5)
ret.bma <- ret.list$obj

## Get the selected genes with probne0 > 0
selected.genes <- ret.list$curr.names[ret.bma$probne0 > 0]

## Get the subset of test data with the genes from the last iteration of bic.surv
curr.test.dat <- testData [, selected.genes]

## Compute the predicted risk scores for the test samples
y.pred.test <- apply (curr.test.dat, 1, predictBicSurv, postprob.vec=ret.bma$postprob, mle.mat=ret.bma$mle)

}
\keyword{survival}



