# Goal:      Carry out iterative bic.surv in R 
# Arguments: x = matrix of independent variables in the training set (variables = columns; samples = rows)
#            surv.time = vector of survival times for all samples in the training set
#	     cens.vec = vector of censor data for all samples in the training set (0 = censored, 1 = uncensored)
#	     curr.mat = matrix of independent variables in the active bic.surv window 
#	     stopVar = 0 to continue iterations, 1 to stop iterations (default 0)
#            nextVar = placeholder indicating the next variable to be brought into the active bic.surv window
#            nbest = number specifying the number of models of each size returned to bic.surv in the BMA package (default 10)
#            maxNvar = size of the active bic.surv window (default 25 variables - may need to be reduced with smaller training sets)
#            maxIter = max number of iterations in repeating bic.surv (default 200000)
#            thresProbne0 = threshold to remove genes with low probne0 (default 1%)
#            verbose = a boolean variable indicating whether or not to print interim information to the console (default FALSE)
# Returns:   updated curr.mat
#            updated nextVar
#            updated stopVar

iterateBMAsurv.train <- function (x, surv.time, cens.vec, curr.mat, stopVar=0, nextVar, nbest=10, maxNvar = 25, maxIter=200000, thresProbne0 = 1, verbose = FALSE, suff.string="") {
	# Iterative bic.surv
	currIter <- 0
	
	while (stopVar == 0 && currIter < maxIter) {
		# Run bic.surv
		if (verbose == TRUE) {
			cat (paste("Current iteration is ", currIter, "\n", sep = ""))
			cat ("Apply bic.surv now\n")
		}
		ret.bic.surv <- bic.surv (x=curr.mat, surv.t=surv.time, cens=cens.vec, nbest=nbest, maxCol=(maxNvar+1))
		
		if (verbose == TRUE) {
			cat ("After bic.surv\n")
		}
		
		# Get a logical vector for which probne0 < thresProbne0; i.e., get a vector of variables whose probne0's are too low
		rmVector <- which (ret.bic.surv$probne0 < thresProbne0)
		if (verbose == TRUE) {
			cat("Posterior Probabilities of selected genes:\n")
			print (ret.bic.surv$probne0)
			cat (paste("Length of rmVector is ", length(rmVector), "\n", sep = ""))
		}
		
		if (any(ret.bic.surv$probne0 < thresProbne0) == FALSE) {
			# No gene to swap in!! Increase threshold
			currMin <- min (ret.bic.surv$probne0)
			if (verbose ==TRUE) {
				cat (paste("No gene to swap! Min probne0 = ", currMin, , "\n", sep=""))
			}
			
			# Assign new threshold
			newThresProbne0 <- currMin + 1
			rmVector <- which (ret.bic.surv$probne0 < newThresProbne0)
			if (verbose == TRUE) {
				cat (paste("New probne0 threshold = ", newThresProbne0, "\n", sep=""))
				cat ("New rmVector after applying increased threshold:\n")
				print (rmVector)
			}
		}

		# Now, any(ret.bic.surv$probne0 < thresProbne0) = TRUE; i.e., there is at least 1 gene to swap, guaranteed
		if (nextVar <= ncol(x)) {
			# Set up new X
			if (verbose == TRUE) {
				cat ("Set up new X\n")
				cat (paste("nextVar is ", nextVar, , "\n", sep=""))
			}
			lastVar <- length(rmVector) + nextVar - 1
						
			# Make sure lastVar <= ncol(x)
			if (lastVar > ncol(x)) {
				rmVector <- rmVector [1: (ncol(x) -nextVar + 1)]
				lastVar <- ncol(x)
			}
			
			# Update curr.mat
			curr.mat[, rmVector] <- x[, nextVar:lastVar]
			
			# Change the colume names as well!!
			dimnames(curr.mat)[[2]][rmVector] <- dimnames(x)[[2]][nextVar:lastVar]
			nextVar <- lastVar + 1
			
		} else {
			# There is no variable to be removed OR exhausted all data
			stopVar <- 1
	        }
	        currIter <- currIter + 1
	}
        
        cat (paste(currIter, ": Explored up to variable # ", nextVar-1, "\n", sep=""))

	# Print out selected genes if iterateBMAsurv has completed
	if (stopVar == 1) {
		cat ("Iterate bic.surv is done!\n")
		cat("Selected genes:\n")
		print (dimnames(curr.mat)[[2]])
		cat("Posterior probabilities of selected genes:\n")
		print (ret.bic.surv$probne0)
		selectedGenes <- which(ret.bic.surv$probne0 >= thresProbne0)
		currfilename <- paste ("final_nbest", nbest, "_adaptThres", thresProbne0, suff.string, ".txt", sep="")
		write.table (curr.mat[, selectedGenes], file=currfilename, sep="\t", quote=FALSE)
	}
	
	list(curr.mat=curr.mat, stopVar=stopVar, nextVar=nextVar)
}

# Goal:      Initialize iterateBMAsurv.train
# Arguments: x = matrix of independent variables in the training set (variables = columns; samples = rows)
#            maxNvar = size of the active bic.surv window (default 25 variables - may need to be reduced with smaller training sets)
# Returns:   curr.mat
#            nextVar
#            stopVar

iterateBMAinit <- function (x, maxNvar = 25) {

	maxNvar <- min (maxNvar, ncol(x))
	curr.mat <- x[, 1:maxNvar]
	stopVar <- 0
	nextVar <- maxNvar + 1

	list(curr.mat=curr.mat, stopVar=stopVar, nextVar=nextVar)
}

# Goal:      Call iterateBMAinit and then iterateBMAsurv.train once
# Caution:   Make sure system has enough memory and that the max number of iterations (maxIter) is high enough
# Arguments: x = matrix of independent variables in the training set (variables = columns; samples = rows) 
#            surv.time = vector of survival times for all samples in the training set
#	     cens.vec = vector of censor data for all samples in the training set (0 = censored, 1 = uncensored)
#            nbest = number specifying the number of models of each size returned to bic.surv in the BMA package (default 10)
#            maxNvar = size of the active bic.surv window (default 25 variables - may need to be reduced with smaller training sets)
#            maxIter = max number of iterations in repeating bic.surv (default 200000)
#            thresProbne0 = threshold to remove genes with low probne0 (default 1%)
#            verbose = a boolean variable indicating whether or not to cat interim information to the console (default FALSE)
# Returns:   A list of 2 components:
# 	     obj = an object of class bic.surv if the algorithm finishes
#            curr.names = dimnames of curr.mat; i.e., input variable names to be passed to bic.surv
#            otherwise, returns -1

iterateBMAsurv.train.wrapper <- function (x, surv.time, cens.vec, nbest=10, maxNvar=25, maxIter=200000, thresProbne0=1, verbose = FALSE, suff.string="") {
	# Get the top "maxNvar" variables
	ret.bma.init <- iterateBMAinit (x, maxNvar)
	
	# Call bic.surv repeatedly
	ret.bma <- iterateBMAsurv.train (x, surv.time, cens.vec, curr.mat=ret.bma.init$curr.mat, stopVar=ret.bma.init$stopVar, nextVar=ret.bma.init$nextVar, nbest, maxNvar, maxIter, thresProbne0, verbose = verbose, suff.string)

	if (ret.bma$stopVar == 1) {
		# Apply bic.surv again using selected genes
		ret.bic.surv <- bic.surv (x=ret.bma$curr.mat, surv.t=surv.time, cens=cens.vec, nbest=nbest, maxCol=(maxNvar+1))
		return (list (obj=ret.bic.surv, curr.names=dimnames(ret.bma$curr.mat)[[2]]))
	} else {
		return (-1)
	}
}

# Goal:      Sort genes in loglik descending order, run iterative bic.surv, predict risk scores in the test set,
#            assign risk categories for test samples, and graph the results
# Arguments: train.dat, test.dat = unsorted training and testing data matrices respectively
# 	     surv.time.train, surv.time.test = survival times for training and test sets respectively
#	     cens.time.train, cens.time.test = censor vectors for training and testing sets respectively (0=censored, 1=uncensored)
#	     p = number of top genes to be used in iterate bic.surv (default 100)
#            nbest = number specifying the number of models of each size returned to bic.surv in the BMA package (default 10)
#            maxNvar = size of the active bic.surv window (default 25 variables - may need to be reduced with smaller training sets)
#            maxIter = max number of iterations in repeating bic.surv (default 200000)
#            thresProbne0 = threshold to remove genes with low probne0 (default 1%)
#            cutPoint = threshold for separating high- from low-risk groups; enter as 50 for a 50% cutoff (default 50)
#            verbose = a boolean variable indicating whether or not to print interim information to the console (default FALSE)
# Returns:   A list with 5 components:
# 	     nvar = number of genes after the final iteration of iterative bic.surv
#	     nmodel = number of models after the final iteration of iterative bic.surv
#	     ypred = predicted risk scores on the test set
#            result.table = tabulated results of the risk groups for the test set
#            statistics = survdiff object containing the outcome statistics

iterateBMAsurv.train.predict.assess <- function  (train.dat, test.dat, surv.time.train, surv.time.test, cens.vec.train, cens.vec.test, p=100, nbest=10, maxNvar=25, maxIter=200000, thresProbne0=1, cutPoint=50, verbose = FALSE, suff.string="") {
 
	# 1. Sort the genes in loglik descending order using Cox Proportional Hazards Regression  
	sorted.genes <- singleGeneCoxph(train.dat, surv.time.train, cens.vec.train)
	sorted.top.genes <- printTopGenes(sorted.genes, p, train.dat)
	currfilename <- paste("sorted_topCoxphGenes_", p, ".txt", sep = "")
	
	# 2. Re-form the training dataset to contain the top p variables in sorted order
	train.dat <- read.table(currfilename, header = TRUE)
	if (verbose == TRUE) {
		cat("Right here, the training data set should be in sorted order:\n")
		print (train.dat[1:min(5, nrow(train.dat)), 1:min(5, ncol(train.dat))])
	}

	# 3. Run iterative bic.surv, starting with calling the iterativeBMAsurv wrapper
	ret.list <- iterateBMAsurv.train.wrapper (x=train.dat, surv.time=surv.time.train, cens.vec=cens.vec.train, nbest, maxNvar, maxIter, thresProbne0, verbose = verbose, suff.string)
	ret.bma <- ret.list$obj

	# 4. Save the results
	write.table (ret.bma$postprob, file=paste("postprob", suff.string, ".txt", sep=""), sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)
	write.table (ret.bma$mle, file=paste("mle", suff.string, ".txt", sep=""), sep="\t", quote=FALSE)
	write.table (ret.bma$probne0, file=paste("probne0", suff.string, ".txt", sep=""), sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)
	write.table (ret.list$curr.names, file=paste("namesx", suff.string, ".txt", sep=""), sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)
	write.table (ret.bma$which, file=paste("which", suff.string, ".txt", sep=""), sep="\t", quote=FALSE)

	# 5. Predict risk scores in the test set
	selected.genes <- ret.list$curr.names [ret.bma$probne0 > 0]
	curr.test.dat <- test.dat [, selected.genes]
	y.pred.test <- apply (curr.test.dat, 1, predictBicSurv, postprob.vec=ret.bma$postprob, mle.mat=ret.bma$mle)

	# 6. Compute risk scores in the training set
	y.pred.train <- apply (train.dat[, selected.genes], 1, predictBicSurv, postprob.vec=ret.bma$postprob, mle.mat=ret.bma$mle)
    	
    	# 7. Assign risk categories for test samples
    	ret.table <- predictiveAssessCategory (y.pred.test, y.pred.train, cens.vec.test, cutPoint)
	cat (paste("# selected genes = ", length(selected.genes), "\n", sep=""))
	cat (paste("# selected models = ", length(ret.bma$postprob), "\n",  sep=""))
	cat ("Risk Table:\n")
	print (ret.table$assign.risk)
	
	risk.groups = ret.table$groups
	if (verbose == TRUE) {
		cat("Risk.groups -- will be error if all samples are assigned to the same risk
		       group or all samples are in the same censor category.\n")
		print (risk.groups)
	}
	
	# 8. Create a survival object from the test set and graph the results
	mySurv.obj <- Surv(surv.time.test, cens.vec.test)
	km.fit <- survfit(mySurv.obj ~ unlist(risk.groups))
	strata.levels <- summary(km.fit)$strata
	filenameThisRun <- paste("km_p_", p, "_nbest_", nbest,"_cutPoint_", cutPoint, "_maxNvar_", maxNvar, ".pdf", sep="")
  
	if ((nrow(ret.table$assign.risk) == 2) && (ncol(ret.table$assign.risk) == 2)) {
	        
	        # There are both high- and low-risk and censored/uncensored samples in the test set; get stats and graph differences
		stats <- survdiff(mySurv.obj ~ unlist(risk.groups))
		print (stats)
		print (stats$var)
		print (stats$chisq)
		pvalue <- 1 - pchisq(stats$chisq, 1)
		
		# Graph	
		pdf(filenameThisRun)
		plot(km.fit, col = c("red", "blue"), main = filenameThisRun, xlab = "Survival Time", ylab = "Survival Probability")
		legend("bottomleft", legend = c(levels(strata.levels)[1], levels(strata.levels)[2]), fill = c("red", "blue"))
		legend("topright", legend = paste("p=", pvalue, sep=""))
		dev.off()
		
		list (nvar=length(selected.genes), nmodel=length(ret.bma$postprob), ypred=y.pred.test, result.table=ret.table$assign.risk, statistics=stats, success=TRUE)
	
	} else {
		# All samples in the test set were assigned as either high-risk or low-risk or
		# all samples were in the same censor category
		cat("ERROR: All test samples were assigned to one risk group or all samples
		       were in the same censor category.\n")
		list(success=FALSE)
	}
	
}

# Goal:      To assess predictive discrimination of bic.surv results by assigning
# 	     samples to 2 discrete risk categories: high and low risk
# Arguments: y.pred.test = predicted risk scores on the test samples
# 	     y.pred.train = computed risk scores on the training samples
#	     cens.vec.test = censored vector on the test set
#            cutPoint = threshold percentage for dividing high- from low-risk (default 50)
# Returns:   the tabulated results of risk categories (row) by censored vector

predictiveAssessCategory <- function (y.pred.test, y.pred.train, cens.vec.test, cutPoint=50) {

	# Define high- and low-risk groups by the empirical argument cutoff point of the risk score from the training set
	# Low-risk group: [min(y.pred.train), ret.quantile[1])
	# High-risk group: [ret.quantile[1], max(y.pred.train)]
	ret.quantile <- quantile (y.pred.train, probs=c(cutPoint/100, 1), na.rm = TRUE)
	
	# Assign each sample to a risk group
	risk.gp.test <- lapply (y.pred.test, assignRiskGroup, quantile.cutPoint=ret.quantile[1])

	# Tabulate results:
	risk.table <- table (as.vector(risk.gp.test, mode="character"), cens.vec.test)
	
	list(assign.risk=risk.table, groups=risk.gp.test)
}

# Goal:      For a given value x, decide which of 2 risk groups x belongs to
# Arguments: x = current sample whose risk group is being assessed
#            quantile.cutPoint = the real-number risk score threshold 
# Returns:   0 for low risk group if x < quantile.cutPoint
#            1 for hi risk group if x >= quantile.cutPoint

assignRiskGroup <- function (x, quantile.cutPoint) {
	if (x < quantile.cutPoint) {
		return ("Low Risk")
	} 
	else {
		return ("High Risk")
	}
}

