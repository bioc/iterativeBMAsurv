\name{iterativeBMAsurv-package}
\alias{iterativeBMAsurv-package}
\alias{iterativeBMAsurv}
\docType{package}
\title{
The Iterative Bayesian Model Averaging (BMA) algorithm for survival analysis
}
\description{
The iterative Bayesian Model Averaging (BMA) algorithm
for survival analysis is a variable selection method
for applying survival analysis to microarray data..
}
\details{
\tabular{ll}{
Package: \tab iterativeBMAsurv\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2008-3-24\cr
License: \tab GPL version 2 or higher\cr
}
The function \code{iterateBMAsurv.train} selects relevant variables by
iteratively applying the \code{bic.surv} function from the \code{BMA} package
until all variables in the training data are exhausted. The variables are 
assumed to be pre-sorted by rank when this function is called. The function 
\code{iterateBMAsurv.train.wrapper} acts as a wrapper for \code{iterateBMAsurv.train}, 
returning the names of the selected variables and an object of class \code{bic.surv} 
if the iterations exhaust all variables in the training set (-1 otherwise). Again, 
the variables are assumed to be pre-sorted by rank, so calling this function 
allows users to experiment with different univariate ranking measures. The function 
\code{iterateBMAsurv.train.predict.assess} combines the training, prediction, and 
assessment phases. It returns a list consisting of the numbers of selected genes 
and models from the training phase, the predicted risk scores of the test samples, 
and the overall survival analysis statistics indicating the difference between risk 
groups (p-value, chi-square statistic, and variance matrix). It also writes a 
Kaplan-Meier survival analysis curve to file, which serves as a pictorial nonparametric
estimator of the difference between risk groups. The variables are not assumed to be 
pre-sorted by rank when this function is called. \code{iterateBMAsurv.train.predict.assess} 
calls \code{singleGeneCoxph}, which ranks the genes based on their log likelihood scores 
using Cox Proportional Hazards Regression. \code{iterateBMAsurv.train.predict.assess} 
calls \code{iterateBMAsurv.train.wrapper} in its training phase, so if Cox Proportional 
Hazards Regression is the desired univariate ranking algorithm, then calling this 
function with the training and testing sets is all that is necessary for a complete 
survival analysis run. The function \code{crossVal} performs k runs of n-fold cross
validation on a training data set, where k and n are specified by the user. 
\code{crossVal} calls \code{iterateBMAsurv.train.predict.assess} during each fold,
so Cox Proportional Hazards Regression is the univariate ranking measure for this
function.
}
\author{
Ka Yee Yeung, University of Washington, Seattle, WA, 
Amalia Annest, University of Washington, Tacoma, WA

Maintainer: Ka Yee Yeung <kayee@u.washington.edu>
             Amalia Annest <amanu@u.washington.edu>
}
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

\seealso{\code{\link{iterateBMAsurv.train}},  
	 \code{\link{iterateBMAsurv.train.wrapper}},
	 \code{\link{iterateBMAsurv.train.predict.assess}},
	 \code{\link{singleGeneCoxph}},
	 \code{\link{predictiveAssessCategory}},
	 \code{\link{crossVal}},
	 \code{\link{trainData}},
	 \code{\link{trainSurv}}, 
	 \code{\link{trainCens}},
	 \code{\link{testData}},
	 \code{\link{testSurv}},
       	 \code{\link{testCens}}
}

\examples{
library (BMA)
library (iterativeBMAsurv)
data(trainData)
data(trainSurv)
data(trainCens)
data(testData)
data(testSurv)
data(testCens)

## Use p=10 genes and nbest=5 for fast computation
ret.bma <- iterateBMAsurv.train.predict.assess (train.dat=trainData, test.dat=testData, surv.time.train=trainSurv, surv.time.test=testSurv, cens.vec.train=trainCens, cens.vec.test=testCens, p=10, nbest=5)

## Extract the statistics from this survival analysis run
number.genes <- ret.bma$nvar
number.models <- ret.bma$nmodel
evaluate.success <- ret.bma$statistics

## Perform 1 run of 2-fold cross validation on the training set, using p=10 genes and nbest=5 for fast computation
cv <- crossVal(exset=trainData, survTime=trainSurv, censor=trainCens, diseaseType="DLBCL", noFolds=2, noRuns=1, p=10, nbest=5)

}
\keyword{multivariate}
\keyword{survival}
