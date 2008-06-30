\name{crossVal}
\alias{crossVal}
\title{Cross Validation for Iterative Bayesian Model Averaging}
\description{This function performs k runs of n-fold cross validation on a 
             training dataset for survival analysis on microarray data, where 
             k and n are specified by the user.}

\usage{crossVal(exset, survTime, censor, diseaseType="cancer", nbest=10, maxNvar=25, p=100, cutPoint=50, verbose=FALSE, noFolds=10, noRuns=10)}

\arguments{
\item{exset}{Data matrix for the training set where columns are variables 
             and rows are observations. In the case of gene expression data, 
             the columns (variables) represent genes, while the rows 
             (observations) represent samples. The data is not assumed
             to be pre-sorted by rank.}
\item{survTime}{Vector of survival times for the patient samples in the 
                training set. Survival times are assumed to be presented 
                in uniform format (e.g., months or days), and the length 
                of this vector should be equal to the number of rows in 
                exset.}
\item{censor}{Vector of censor data for the patient samples in the 
              training set. In general, 0 = censored and 1 = uncensored. 
              The length of this vector should equal the number of rows 
              in exset and the number of elements in survTime.}
\item{diseaseType}{String denoting the type of disease in the training dataset 
                   (used for writing to file). Default is 'cancer'.}
\item{nbest}{A number specifying the number of models of each size 
             returned to \code{bic.surv} in the \code{BMA} package. 
	     The default is 10.}
\item{maxNvar}{A number indicating the maximum number of variables used in
	       each iteration of \code{bic.surv} from the \code{BMA} package.
	       The default is 25.}
\item{p}{A number indicating the maximum number of top univariate genes
	 used in the iterative \code{bic.surv} algorithm.  This number is 
	 assumed to be less than the total number of genes in the training data.
	 A larger p usually requires longer computational time as more iterations 
	 of the \code{bic.surv} algorithm are potentially applied. The default is 
	 100.}
\item{cutPoint}{Threshold percent for separating high- from low-risk groups. 
                The default is 50.}
\item{verbose}{A boolean variable indicating whether or not to print interim 
               information to the console. The default is FALSE.}
\item{noFolds}{A number specifying the desired number of folds in each cross
               validation run. The default is 10.}
\item{noRuns}{A number specifying the desired number of cross validation runs.
              The default is 10.}		    
}

\details{This function performs k runs of n-fold cross validation, where k and n
         are specified by the user through the \code{noRuns} and \code{noFolds} 
         arguments respectively. For each run of cross validation, the training
         set, survival times, and censor data are re-ordered according to a 
         random permutation. For each fold of cross validation, \eqn{1/n}th of the data 
         is set aside to act as the validation set. In each fold, the 
         \code{iterateBMAsurv.train.predict.assess} function is called in order 
         to carry out a complete run of survival analysis. This means the 
         univariate ranking measure for this cross validation function is Cox 
         Proportional Hazards Regression; see \code{iterateBMAsurv.train.wrapper} 
         to experiment with alternate univariate ranking methods. With each run 
         of cross validation, the survival analysis statistics are saved and 
         written to file.}

\value{The output of this function is a series of files giving information on
       cross validation results. The file beginning with 'foldresults' contains
       information for every fold in the form of a 2 x 2 table indicating the 
       number of test samples in each category (high-risk or censored, 
       high-risk or uncensored, low-risk or censored, low-risk or uncensored). This file
       also gives the accumulated percentage of  uncensored statistic from each run. The file
       beginning with 'runresults' gives the total number of test samples assigned
       to each category along with percentage uncensored across the entire run. The end of 
       this file contains this same information, averaged across all runs. The 
       file beginning with 'stats' gives the statistics from each fold, including 
       the p-value, chi-square statistic, and variance matrix. Finally, the file 
       beginning with 'avg\_p\_value\_chi\_square' gives the overall means and standard 
       deviations of the p-values and chi-square statistics across all runs and all 
       folds.}


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
\note{The \code{BMA} package is required. Also, smaller training sets may lead to
      cross validation folds where all test samples are assigned to one risk group 
      or all samples are in the same censor category (all samples are either censored
      or uncensored). In this case, the fold is skipped, and cross validation proceeds 
      from the next fold. This particular error will be evidenced by a missing fold 
      result in the output files. All averages will be calculated as if this fold had 
      never occurred.}

\seealso{\code{\link{iterateBMAsurv.train.predict.assess}}
	 \code{\link{iterateBMAsurv.train.wrapper}},  
	 \code{\link{iterateBMAsurv.train}},
	 \code{\link{singleGeneCoxph}},
	 \code{\link{predictBicSurv}},
	 \code{\link{predictiveAssessCategory}},
	 \code{\link{trainData}},
	 \code{\link{trainSurv}}, 
	 \code{\link{trainCens}}
	 }

\examples{
library (BMA)
library(iterativeBMAsurv)
data(trainData)
data(trainSurv)
data(trainCens)

## Perform 1 run of 2-fold cross validation on the training set, using p=10 genes and nbest=5 for fast computation
cv <- crossVal (exset=trainData, survTime=trainSurv, censor=trainCens, diseaseType="DLBCL", noRuns=1, noFolds=2, p=10, nbest=5)

## Upon completion of this function, all relevant output files will be in the working R directory.

}
\keyword{survival}



