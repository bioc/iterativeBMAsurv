\name{iterateBMAsurv.train.predict.assess}
\alias{iterateBMAsurv.train.predict.assess}
\title{Iterative Bayesian Model Averaging: training, prediction, assessment}
\description{Survival analysis and variable selection on microarray data.
	     This is a multivariate technique to select a small number
	     of relevant variables (typically genes) to perform survival 
	     analysis on microarray data. This function performs the training,
	     prediction, and assessment steps.  The data is not assumed to be
	     pre-sorted by rank before this function is called.}

\usage{iterateBMAsurv.train.predict.assess (train.dat, test.dat, surv.time.train, surv.time.test, cens.vec.train, cens.vec.test, p=100, nbest=10, maxNvar=25, maxIter=200000, thresProbne0=1, cutPoint=50, verbose = FALSE, suff.string="")}

\arguments{
\item{train.dat}{Data matrix for the training set where columns are variables 
                 and rows are observations. In the case of gene expression data, 
                 the columns (variables) represent genes, while the rows 
                 (observations) represent samples. If Cox Proportional Hazards 
                 Regression is the desired univariate ranking measure, the data 
                 does not need to be pre-sorted. To use a different univariate 
                 ranking measure, see \code{iterateBMAsurv.train.wrapper}.}
\item{test.dat}{Data matrix for the test set where columns are variables 
                and rows are observations. In the case of gene expression data, 
                the columns (variables) represent genes, while the rows 
                (observations) represent samples. The test set should contain
                the same variables as the training set.}
\item{surv.time.train}{Vector of survival times for the patient samples in the 
                       training set. Survival times are assumed to be presented 
                       in uniform format (e.g., months or days), and the length 
                       of this vector should be equal to the number of rows in 
                       train.dat.}
\item{surv.time.test}{Vector of survival times for the patient samples in the 
                      test set. Survival times are assumed to be presented 
                      in uniform format (e.g., months or days), and the length 
                      of this vector should be equal to the number of rows in 
                      test.dat.}
\item{cens.vec.train}{Vector of censor data for the patient samples in the 
                      training set. In general, 0 = censored and 1 = uncensored. 
                      The length of this vector should equal the number of rows 
                      in train.dat and the number of elements in surv.time.train.}
\item{cens.vec.test}{Vector of censor data for the patient samples in the 
                     test set. In general, 0 = censored and 1 = uncensored. 
                     The length of this vector should equal the number of rows 
                     in train.dat and the number of elements in surv.time.test.}
\item{p}{A number indicating the maximum number of top univariate genes
	 used in the iterative \code{bic.surv} algorithm.  This number is 
	 assumed to be less than the total number of genes in the training data.
	 A larger p usually requires longer computational time as more iterations 
	 of the \code{bic.surv} algorithm are potentially applied. The default is 
	 100.}
\item{nbest}{A number specifying the number of models of each size 
             returned to \code{bic.surv} in the \code{BMA} package. 
	     The default is 10.}
\item{maxNvar}{A number indicating the maximum number of variables used in
	       each iteration of \code{bic.surv} from the \code{BMA} package.
	       The default is 25.}
\item{maxIter}{A number indicating the maximum iterations of \code{bic.surv}. 
               The default is 200000.}
\item{thresProbne0}{A number specifying the threshold for the posterior
                    probability that each variable (gene) is non-zero (in
		    percent).  Variables (genes) with such posterior 
		    probability less than this threshold are dropped in
		    the iterative application of \code{bic.surv}. The default
		    is 1 percent.}
\item{cutPoint}{Threshold percent for separating high- from low-risk groups. 
                The default is 50.}
\item{verbose}{A boolean variable indicating whether or not to print interim 
               information to the console. The default is FALSE.}
\item{suff.string}{A string for writing to file.}
}

\details{This function consists of the training phase, the prediction phase,
         and the assessment phase. The training phase first orders all the 
         variables (genes) by a univariate measure called Cox Proportional 
         Hazards Regression, and then iteratively applies the \code{bic.surv} 
         algorithm from the \code{BMA} package. The prediction phase uses the 
         variables (genes) selected in the training phase to predict the risk 
         scores of the patient samples in the test set. In the assessment phase, 
         the risk scores are used to designate each test sample as either 
         high-risk or low-risk based on the user-designated \code{cutPoint}. 
         Prediction accuracy is measured by the p-value difference between 
         groups as calculated through the central chi-square distribution. In
         addition, a Kaplan-Meier Survival Analysis Curve illustrating the
         difference between risk groups is written to file in the working R
         directory. If Cox Proportional Hazards Regression is the desired 
         univariate ranking algorithm, then calling this function with the 
         training and testing sets is all that is necessary for a complete 
         survival analysis run.}

\value{If all samples are assigned to a single risk group or all samples are in
       the same censor category, an error message is printed and a boolean variable 
       \code{success} is returned as FALSE.If both risk groups are present in the 
       patient test samples, a Kaplan-Meier Survival Analysis Curve is written to file, 
       and a list with 6 components is returned:
\item{nvar}{The number of variables selected by the last iteration of \code{bic.surv}.}
\item{nmodel}{The number of models selected by the last iteration of \code{bic.surv}.}
\item{ypred}{The predicted risk scores on the test samples.}
\item{result.table}{A 2 x 2 table indicating the number of test samples in each 
                    category (high-risk/censored, high-risk/uncensored, 
                    low-risk/censored, low-risk/uncensored).}
\item{statistics}{An object of class \code{survdiff} that contains the statistics
                  from survival analysis, including the variance matrix, chi-square
                  statistic, and p-value.}
\item{success}{A boolean variable returned as TRUE if both risk groups are present
               in the patient test samples.}
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
\note{The \code{BMA} package is required.}

\seealso{\code{\link{iterateBMAsurv.train.wrapper}},  
	 \code{\link{iterateBMAsurv.train}},
	 \code{\link{singleGeneCoxph}},
	 \code{\link{predictBicSurv}},
	 \code{\link{predictiveAssessCategory}},
	 \code{\link{trainData}},
	 \code{\link{trainSurv}}, 
	 \code{\link{trainCens}},
         \code{\link{testData}},
         \code{\link{testSurv}},
       	 \code{\link{testCens}}
}

\examples{
library (BMA)
library(iterativeBMAsurv)
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

}
\keyword{multivariate}
\keyword{survival}



