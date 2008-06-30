\name{iterateBMAsurv.train.wrapper}
\alias{iterateBMAsurv.train.wrapper}
\title{Iterative Bayesian Model Averaging: training}
\description{This function is a wrapper for \code{iterateBMAsurv.train}, which 
             repeatedly calls \code{bic.surv} from the \code{BMA} package 
             until all variables are exhausted. At the point when this function
             is called, the variables in the dataset are assumed to be 
             pre-sorted by rank.}
	     
\usage{iterateBMAsurv.train.wrapper (x, surv.time, cens.vec, nbest=10, maxNvar=25, maxIter=200000, thresProbne0=1, verbose=FALSE, suff.string="")}

\arguments{
\item{x}{Data matrix where columns are variables and rows are observations.  
         The variables (columns) are assumed to be sorted using a univariate 
         measure. In the case of gene expression data, the columns (variables) 
	 represent genes, while the rows (observations) represent samples.}
\item{surv.time}{Vector of survival times for the patient samples. Survival times 
                 are assumed to be presented in uniform format (e.g., months or 
                 days), and the length of this vector should be equal to the number 
                 of rows in x.}
\item{cens.vec}{Vector of censor data for the patient samples. In general,
                0 = censored and 1 = uncensored. The length of this vector
                should equal the number of rows in x and the number of elements 
                in surv.time.}
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
\item{verbose}{A boolean variable indicating whether or not to print interim 
               information to the console. The default is FALSE.}
\item{suff.string}{A string for writing to file.}
}

\details{In this wrapper function for \code{iterateBMAsurv.train}, the variables 
         are assumed to be sorted, and \code{bic.surv} is called repeatedly 
         until all the variables have been exhausted.  In the first application 
         of the \code{bic.surv} algorithm, the top {\tt maxNvar} univariate 
         ranked genes are used. After each application of the \code{bic.surv} 
         algorithm, the genes with {\tt probne0} < {\tt thresProbne0}
	 are dropped, and the next univariate ordered genes are added
	 to the \code{bic.surv} window. The function 
	 \code{iterateBMAsurv.train.predict.assess} calls \code{SingleGeneCoxph} 
	 before calling this function. Using this function directly, users can 
	 experiment with alternative univariate measures.}

\value{If {\tt maxIter} is reached or the iterations stop before all variables 
       are exhausted, -1 is returned. If all variables are exhausted, two items 
       are returned: 
\item{curr.names}{A vector containing the names of the variables (genes)
                  from the last iteration of \code{bic.surv}}. 
\item{obj}{An object of class \code{bic.surv} returned by the last iteration of 
           \code{bic.surv}. The object of class \code{bic.surv} is a list 
           consisting of the following components:
        \item{namesx}{the names of the variables in the last iteration of 
                      \code{bic.surv}.}
	\item{postprob}{The posterior probabilities of the models selected.}
	\item{label}{Labels identifying the models selected.}
	\item{bic}{Values of BIC for the models.}
	\item{size}{The number of independent variables in each of the models.}
	\item{which}{A logical matrix with one row per model and one column per 
		     variable indicating whether that variable is in the model.}
	\item{probne0}{The posterior probability that each variable is non-zero 
		       (in percent).}
	\item{postmean}{The posterior mean of each coefficient (from model averaging).}
	\item{postsd}{The posterior standard deviation of each coefficient 
		      (from model averaging).}
	\item{condpostmean}{The posterior mean of each coefficient conditional on 
			    the variable being included in the model.}
	\item{condpostsd}{The posterior standard deviation of each coefficient 
			  conditional on the variable being included in the model.}
	\item{mle}{Matrix with one row per model and one column per variable giving 
		   the maximum likelihood estimate of each coefficient for each model.}
	\item{se}{Matrix with one row per model and one column per variable giving 
		  the standard error of each coefficient for each model.}
	\item{reduced}{A logical indicating whether any variables were dropped 
		       before model averaging.}
	\item{dropped}{A vector containing the names of those variables dropped 
		       before model averaging.}
	\item{call}{The matched call that created the bma.lm object.}}
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

\seealso{\code{\link{iterateBMAsurv.train.predict.assess}},  
	 \code{\link{iterateBMAsurv.train}},
	 \code{\link{predictiveAssessCategory}},
	 \code{\link{singleGeneCoxph}},
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

## Training data should be pre-sorted before beginning

## Run iterative bic.surv, using nbest=5 for fast computation
ret.list <- iterateBMAsurv.train.wrapper (x=trainData, surv.time=trainSurv, cens.vec=trainCens, nbest=5)

## Extract the \code{bic.surv} object
ret.bma <- ret.list$obj

## Extract the names of the genes from the last iteration of \code{bic.surv}
gene.names <- ret.list$curr.names

}
\keyword{multivariate}
\keyword{survival}



