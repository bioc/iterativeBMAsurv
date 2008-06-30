\name{printTopGenes}
\alias{printTopGenes}
\title{Write a training set including the top-ranked G variables from a sorted matrix to file}
\description{This function takes a matrix of rank-ordered variables and writes a training set
             containing the top G variables in the matrix to file.}
             
\usage{printTopGenes (retMatrix, numGlist=c(10, 30, 50, 100, 500, 1000, ncol(trainData)), trainData, myPrefix="sorted_topCoxphGenes_")}

\arguments{
\item{retMatrix}{A three-column matrix where the first column contains the sorted variable 
                 names (the top log-ranked variable appears first), the second column contains
                 the original index of the variables, and the third column contains the
                 variable ranking from 1 to ncol(trainData).}
\item{numGlist}{A list of values for the desired number of top-ranked variables to be
                written to file. A separate file will be written for each number 
                G in the list, containing genes 1:G (default = c(10, 30, 50, 
                100, 500, 1000, ncol(trainData))).}
\item{trainData}{Data matrix where columns are variables and rows are observations.
                 In the case of gene expression data, the columns (variables) 
                 represent genes, while the rows (observations) represent patient 
                 samples.}
\item{myPrefix}{A string prefix for the filename (default = 'sorted\_topCoxphGenes\_').}
}

\details{This function is called by \code{iterateBMAsurv.train.predict.assess}. It is meant
         to be used in conjunction with \code{singleGeneCoxph}, as the \code{retMatrix} 
         argument is returned by \code{singleGeneCoxph}.}

\value{A file or files consisting of the training data sorted in descending order
       by the top-ranked G variables (one file for each G in numGList).}


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
         \code{\link{singleGeneCoxph}},
         \code{\link{trainData}},
         \code{\link{trainSurv}}, 
         \code{\link{trainCens}}, 
}

\examples{
library(BMA)
library(iterativeBMAsurv)
data(trainData)
data(trainSurv)
data(trainCens)

## Start by ranking and sorting the genes; in this case we use the Cox Proportional Hazards Model
sorted.genes <- singleGeneCoxph(trainData, trainSurv, trainCens)

## Write top 100 genes to file
sorted.top.genes <- printTopGenes(retMatrix=sorted.genes, 100, trainData)

## The file, 'sorted_topCoxphGenes_100', is now in the working R directory.

}
\keyword{univar}
\keyword{print}


