\name{iterateBMAsurv.train}
\alias{iterateBMAsurv.train}
\title{Iterative Bayesian Model Averaging: training}
\description{Survival analysis and variable selection on microarray data.
	     This is a multivariate technique to select a small number
	     of relevant variables (typically genes) to perform survival 
	     analysis on microarray data.  This function performs the 
	     training phase. It repeatedly calls \code{bic.surv} from the
	     \code{BMA} package until all variables are exhausted. The 
	     variables in the dataset are assumed to be pre-sorted by rank.}
	     
\usage{iterateBMAsurv.train (x, surv.time, cens.vec, curr.mat, stopVar=0, nextVar, nbest=10, maxNvar=25, maxIter=200000, thresProbne0=1, verbose = FALSE, suff.string="")}

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
\item{curr.mat}{Matrix of independent variables in the active \code{bic.surv}
                window. There can be at most \code{maxNvar} variables in the 
                window at any given time.}
\item{stopVar}{0 to continue iterations, 1 to stop iterations (default 0)}
\item{nextVar}{Integer placeholder indicating the next variable to be brought 
               into the active \code{bic.surv} window.}
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

\details{The training phase consists of first ordering all the variables 
         (genes) by a univariate measure such as Cox Proportional Hazards
         Regression, and then iteratively applying the \code{bic.surv} algorithm
	 from the \code{BMA} package.  In the first application of
	 the \code{bic.surv} algorithm, the top \code{maxNvar} univariate
	 ranked genes are used.  After each application of the \code{bic.surv}
	 algorithm, the genes with \code{probne0} < \code{thresProbne0}
	 are dropped, and the next univariate ordered genes are added
	 to the active \code{bic.surv} window.}

\value{On the last iteration of \code{bic.surv}, four items are returned:
\item{curr.mat}{A vector containing the names of the variables (genes) 
                from the final iteration of \code{bic.surv}}.
\item{stopVar}{The ending value of stopVar after all iterations.}
\item{nextVar}{The ending value of nextVar after all iterations.}
\item{}{An object of class \code{bic.surv} resulting from the last 
           iteration of \code{bic.surv}.  The object is a list consisting 
           of the following components:
    \describe{
        \item{namesx}{the names of the variables in the last iteration of 
                      \code{bic.surv}.}
	\item{postprob}{the posterior probabilities of the models selected.}
	\item{label}{labels identifying the models selected.}
	\item{bic}{values of BIC for the models.}
	\item{size}{the number of independent variables in each of the models.}
	\item{which}{a logical matrix with one row per model and one column per 
		     variable indicating whether that variable is in the model.}
	\item{probne0}{the posterior probability that each variable is non-zero 
		       (in percent).}
	\item{postmean}{the posterior mean of each coefficient (from model averaging).}
	\item{postsd}{the posterior standard deviation of each coefficient 
		      (from model averaging).}
	\item{condpostmean}{the posterior mean of each coefficient conditional on 
			    the variable being included in the model.}
	\item{condpostsd}{the posterior standard deviation of each coefficient 
			  conditional on the variable being included in the model.}
	\item{mle}{matrix with one row per model and one column per variable giving 
		   the maximum likelihood estimate of each coefficient for each model.}
	\item{se}{matrix with one row per model and one column per variable giving 
		  the standard error of each coefficient for each model.}
	\item{reduced}{a logical indicating whether any variables were dropped 
		       before model averaging.}
	\item{dropped}{a vector containing the names of those variables dropped 
		       before model averaging.}
	\item{call}{the matched call that created the bma.lm object.}
     }
  }
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
	 \code{\link{iterateBMAsurv.train.predict.assess}},
	 \code{\link{singleGeneCoxph}},
	 \code{\link{predictBicSurv}},
	 \code{\link{trainData}},
	 \code{\link{trainSurv}}, 
         \code{\link{trainCens}},
         \code{\link{testData}}
         }

\examples{
library(BMA)
library(iterativeBMAsurv)
data(trainData)
data(trainSurv)
data(trainCens)
data(testData)

## Training data should be pre-sorted before beginning

## Initialize the matrix for the active bic.surv window with variables 1 through maxNvar
maxNvar <- 25
curr.mat <- trainData[, 1:maxNvar]
nextVar <- maxNvar + 1

## Training phase: select relevant genes using nbest=5 for fast computation
ret.bic.surv <- iterateBMAsurv.train (x=trainData, surv.time=trainSurv, cens.vec=trainCens, curr.mat, stopVar=0, nextVar, nbest=5, maxNvar=25)

# Apply bic.surv again using selected genes
ret.bma <- bic.surv (x=ret.bic.surv$curr.mat, surv.t=trainSurv, cens=trainCens, nbest=5, maxCol=(maxNvar+1))

## Get the matrix for genes with probne0 > 0
ret.gene.mat <- ret.bic.surv$curr.mat[ret.bma$probne0 > 0]

## Get the gene names from ret.gene.mat
selected.genes <- dimnames(ret.gene.mat)[[2]]

## Show the posterior probabilities of selected models
ret.bma$postprob

## Get the subset of test data with the genes from the last iteration of
## 'bic.surv'
curr.test.dat <- testData[, selected.genes]

## Compute the predicted risk scores for the test samples
y.pred.test <- apply (curr.test.dat, 1, predictBicSurv, postprob.vec=ret.bma$postprob, mle.mat=ret.bma$mle)
}

\keyword{multivariate}
\keyword{survival}



