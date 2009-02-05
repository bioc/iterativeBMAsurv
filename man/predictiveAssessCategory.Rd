\name{predictiveAssessCategory}
\alias{predictiveAssessCategory}
\title{Risk Groups: assignment of patient test samples}
\description{This function assigns a risk group (high-risk or low-risk) to each
             patient sample in the test set based on the value of the patient's
             predicted risk score. The \code{cutPoint} between high-risk and 
             low-risk is designated by the user.}
             
\usage{predictiveAssessCategory (y.pred.test, y.pred.train, cens.vec.test, cutPoint=50)}

\arguments{
\item{y.pred.test}{A vector containing the predicted risk scores of the test samples.}
\item{y.pred.train}{A vector containing the computed risk scores of the training samples.}
\item{cens.vec.test}{A vector of censor data for the patient samples in the 
                     test set. In general, 0 = censored and 1 = uncensored.}
\item{cutPoint}{Threshold percent for separating high- from low-risk groups. 
                The default is 50.}	                     
}

\details{This function begins by using the computed risk scores of the training
         set (\code{y.pred.train}) to define a real-number empirical cutoff point
         between high- and low-risk groups. The cutoff point is determined by
         the percentile \code{cutPoint} as designated by the user. The predicted
         risk scores from the test samples are then matched against this cutoff
         point to determine whether they belong in the high-risk or the low-risk
         category.}

\value{A list consisting of 2 components:
\item{assign.risk}{A 2 x 2 table indicating the number of test samples in each 
                   category (high-risk/censored, high-risk/uncensored, 
                   low-risk/censored, low-risk/uncensored).}
\item{groups}{A list of all patient samples in the test set with their 
              corresponding 'High-risk' or 'Low-risk' designations.}
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
\note{}

\seealso{\code{\link{iterateBMAsurv.train.predict.assess}},
         \code{\link{predictBicSurv}},
         \code{\link{singleGeneCoxph}},
         \code{\link{printTopGenes}},
         \code{\link{trainData}},
	 \code{\link{trainSurv}}, 
	 \code{\link{trainCens}},
         \code{\link{testData}},
         \code{\link{testSurv}}, 
	 \code{\link{testCens}},
}

\examples{
library(BMA)
library(iterativeBMAsurv)
data(trainData)
data(trainSurv)
data(trainCens)
data(testData)
data(testSurv)
data(testCens)

## Training should be pre-sorted before beginning

## Initialize the matrix for the active bic.surv window with variables 1 through maxNvar
maxNvar <- 25
curr.mat <- trainData[, 1:maxNvar]
nextVar <- maxNvar + 1

## Training phase: select relevant genes, using nbest=5 for fast computation
ret.bic.surv <- iterateBMAsurv.train (x=trainData, surv.time=trainSurv, cens.vec=trainCens, curr.mat, stopVar=0, nextVar, maxNvar=25, nbest=5)

# Apply bic.surv again using selected genes
ret.bma <- bic.surv (x=ret.bic.surv$curr.mat, surv.t=trainSurv, cens=trainCens, nbest=5, maxCol=(maxNvar+1))

## Get the matrix for genes with probne0 > 0
ret.gene.mat <- ret.bic.surv$curr.mat[ret.bma$probne0 > 0]

## Get the gene names from ret.gene.mat
selected.genes <- dimnames(ret.gene.mat)[[2]]

## Show the posterior probabilities of selected models
ret.bma$postprob

## Get the subset of test data with the genes from the last iteration of 'bic.surv'
curr.test.dat <- testData[, selected.genes]

## Compute the predicted risk scores for the test samples
y.pred.test <- apply (curr.test.dat, 1, predictBicSurv, postprob.vec=ret.bma$postprob, mle.mat=ret.bma$mle)

## Compute the risk scores in the training set
y.pred.train <- apply (trainData[, selected.genes], 1, predictBicSurv, postprob.vec=ret.bma$postprob, mle.mat=ret.bma$mle)

## Assign risk categories for test samples
ret.table <- predictiveAssessCategory (y.pred.test, y.pred.train, testCens, cutPoint=50) 

## Extract risk group vector and risk group table
risk.list <- ret.table$groups
risk.table <- ret.table$assign.risk

## Create a survival object from the test set
mySurv.obj <- Surv(testSurv, testCens)

## Extract statistics including p-value and chi-square
stats <- survdiff(mySurv.obj ~ unlist(risk.list))

## The entire block of code above can be executed simply by calling
## 'iterateBMAsurv.train.predict.assess' 

}
\keyword{survival}




