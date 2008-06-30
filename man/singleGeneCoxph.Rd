\name{singleGeneCoxph}
\alias{singleGeneCoxph}
\title{Univariate Cox Proportional Hazards Model for selecting top log-ranked predicitve variables}
\description{This is a univariate technique to rank variables by their predictive 
             relevance for use in survival analysis on microarray data. The 
             log likelihood is computed for each indiviual variable, where a larger 
             log likelihood value indicates a higher rank.}
             
\usage{singleGeneCoxph(trainData, survData, censoredData)}

\arguments{
\item{trainData}{Data matrix where columns are variables and rows are observations.
                 In the case of gene expression data, the columns (variables) 
                 represent genes, while the rows (observations) represent patient 
                 samples.}
\item{survData}{Vector of survival times for the patient samples. Survival times 
                are assumed to be presented in uniform format (e.g., months or 
                days), and the length of this vector should be equal to the number 
                of rows in trainData.}
\item{censoredData}{Vector of censor data for the patient samples. In general,
                    0 = censored and 1 = uncensored. The length of this vector
                    should equal the number of rows in trainData and the number
                    of elements in survData.}
}

\details{This function is called by \code{iterateBMAsurv.train.predict.assess}.}

\value{This function returns a sorted three-column matrix of the training data variables.
       The first column gives the variable names with the top log-ranked variable appearing
       first. The second column gives the original indexes of the variables, and the third
       column gives the rank of the variables from 1 through ncol(trainData). The matrix is 
       also written to file in the working R directory under the filename 'sorted\_loglik.txt'.}


\references{
Annest, A., Yeung, K.Y., Bumgarner, R.E., and Raftery, A.E. (2008).
Iterative Bayesian Model Averaging for Survival Analysis.
Manuscript in Progress.

Cox, D. (1972). Regression Models and Life Tables. 
Journal of the Royal Statistical Society Series B 34: 187-220.

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
         \code{\link{printTopGenes}},
         \code{\link{trainData}},
         \code{\link{trainSurv}}, 
         \code{\link{trainCens}}
}

\examples{
library(BMA)
library(iterativeBMAsurv)
data(trainData)
data(trainSurv)
data(trainCens)

sorted.genes <- singleGeneCoxph(trainData, trainSurv, trainCens)

## Write top 100 genes to file
sorted.top.genes <- printTopGenes(retMatrix=sorted.genes, 100, trainData)

## The file, 'sorted_topCoxphGenes_100', is now in the working R directory.

}
\keyword{univar}
\keyword{regression}


