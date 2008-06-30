# Goal:      Complete k runs of n-fold cross validation on a training dataset
# Requires   BMA, iterativeBMA, and iterativeBMAsurv libraries
# Arguments: exset = matrix of independent variables in the training set (variables = columns; samples = rows)
#            survTime = vector of survival times for all samples in the training set
#            censor = vector of censor data for all samples in the training set (0 = censored, 1 = uncensored)
#            diseaseType = string denoting the type of disease in the training dataset (used for writing to file; default "cancer")
#	     nbest = number specifying the number of models of each size returned to bic.surv in the BMA package (default 10)
#            maxNvar = size of the active bic.surv window (default 25 - may need to be reduced with smaller training sets)
#            p = number of top-ranked genes to be included in the iterations (default 100)
#            cutPoint = threshold for separating high- from low-risk groups; enter as 50 for a 50% cutoff (default 50)
#            verbose = a boolean variable indicating whether or not to print  interim information to the console (default FALSE)
#            noFolds = the number of folds in each cross validation run (default 10)
#            noRuns = the number of cross validation runs (default 10)
# Output:    files giving fold results, run results, stats-by-fold results, and average-over-all-runs results
             
crossVal <- function(exset, survTime, censor, diseaseType="cancer", nbest=10, maxNvar=25, p=100, cutPoint=50, verbose = FALSE, noFolds=10, noRuns=10) {

	# Initialize matrices and filenames
	cv <- crossVal.init(exset, noFolds, noRuns, nbest, maxNvar, p, cutPoint)
	
	# Perform cross validation runs from 1:noRuns
	for (i in 1:noRuns) {
		cv.run <- crossVal.run(exset, survTime, censor, diseaseType, nbest, maxNvar, p, cutPoint, cv$rp.mat, cv$overall.average.results, cv$p.val.result.mat, cv$chi.sq.result.mat, noFolds, cv, run=i, verbose=verbose)
		cv$overall.average.results <- cv.run$overall.average.results
		cv$rp.mat <- cv.run$rp.mat
		cv$p.val.result.mat <- cv.run$p.val.result.mat
		cv$chi.sq.result.mat <- cv.run$chi.sq.result.mat
	}
		
	# Calculate overall average result matrix and final uncensored percentage
	calc <- crossVal.final.calc(cv$overall.average.results, noRuns)
	cv$overall.average.results <- calc$overall.average.results
	
	# Make sure there are no zero-values in the p-value and chi-square matrices for mean/sd calculations
	non.zero.p.val.mat <- cv$p.val.result.mat[cv$p.val.result.mat > 0]
	non.zero.chi.square.mat <- cv$chi.sq.result.mat[cv$chi.sq.result.mat > 0]
	
	# Print and write final overall CV results to file
	cat("Overall average result matrix over all runs:\n")
	print (cv$overall.average.results)
	write.table("Overall average results across all CV runs:", cv$filenamerun, append=TRUE, quote=FALSE)
	write.table(cv$overall.average.results, cv$filenamerun, append=TRUE, quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
	
	# Print and write final overall p-value and chi-square statistics to file
	cat("Average p-value across all folds and all runs:\n")
	print (mean(non.zero.p.val.mat))
	cat("Standard deviation of p-values across all folds and all runs:\n")
	print (sd(non.zero.p.val.mat))
	cat("Average chi-square value across all folds and all runs:\n")
	print (mean(cv$chi.sq.result.mat))
	cat("Standard deviation of chi-square across all folds and all runs:\n")
	print (sd(non.zero.chi.square.mat))
	statsfilename <- paste("avg_p_value_chi_square_", p, "_genes_nbest_", nbest, "_maxNvar_", maxNvar, "_cutPoint_", cutPoint, ".txt", sep="")  
	write.table("P-value matrix:", statsfilename, quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE)
	write.table(cv$p.val.result.mat, statsfilename, quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE)
	write.table("Mean:", statsfilename, quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE)
	write.table(mean(non.zero.p.val.mat), statsfilename, quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE)
	write.table("Standard Deviation:", statsfilename, quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE)
	write.table(sd(non.zero.p.val.mat), statsfilename, quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE)
	write.table("Chi-Square Matrix:", statsfilename, quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE)
	write.table(cv$chi.sq.result.mat, statsfilename, quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE)
	write.table("Mean:", statsfilename, quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE)
	write.table(mean(non.zero.chi.square.mat), statsfilename, quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE)
	write.table("Standard Deviation:", statsfilename, quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE)
	write.table(sd(non.zero.chi.square.mat), statsfilename, quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE)
}

# Goal:      Initialize matrices and file names for cross validation runs
# Arguments: exset = matrix of independent variables in the training set (variables = columns; samples = rows) 
#            noFolds = the number of folds in each cross validation run (default 10)
#            noRuns = the number of cross validation runs (default 10)
#            nbest = number specifying the number of models of each size returned to bic.surv in the BMA package (default 10)
#            maxNvar = size of the active bic.surv window (default 25 - may need to be reduced with smaller training sets)
#            p = number of top-ranked genes to be included in the iterations (default 100)
#            cutPoint = threshold for separating high- from low-risk groups; enter as 50 for a 50% cutoff (default 50)
# Returns:   noSamples = number of samples in the training set
#            noGenes = number of genes in the training set
#            rp.mat = matrix to store the permutation order for each run
#            overall.average.results = matrix to store the overall high-risk/low-risk life/death information
#            p.val.result.mat = matrix to store p-value results for all runs
#            chi_sqaure_result_matrix = matrix to store chi square results for all runs
#            filenamerun = filename for the run results
#            filenamefold = filename for the fold results

crossVal.init <- function(exset, noFolds, noRuns, nbest, maxNvar, p, cutPoint) {
	
	# Save the number of patient samples and the number of genes
	noSamples <- nrow(exset)
	noGenes <- ncol(exset)

	# Random permutation matrix - CV runs are rows and samples are columns
	rp.mat <- mat.or.vec(noRuns, noSamples)

	# Average results matrix to average over all runs
	overall.average.results <- mat.or.vec(2, 3)
	rownames(overall.average.results) <- c("Low Risk", "High Risk")
	colnames(overall.average.results) <- c("Censored", "Uncensored", "Percent Uncensored")

	# Matrices to store p-values and chi-square statistics - CV runs are rows and folds are columns
	p.val.result.mat <- mat.or.vec(noRuns, noFolds)
	chi.sq.result.mat <- mat.or.vec(noRuns, noFolds)

	# Create filenames to store run- and fold-results
	filenamerun = paste("runresults_", p, "_genes_nbest_", nbest, "_maxNvar_", maxNvar, "_cutPoint_", cutPoint, ".txt", sep="")  
	filenamefold = paste("foldresults_", p, "_genes_nbest_", nbest, "_maxNvar_", maxNvar, "_cutPoint_", cutPoint, ".txt", sep="")  
	
	list(noSamples=noSamples, noGenes=noGenes, rp.mat=rp.mat, overall.average.results=overall.average.results, p.val.result.mat=p.val.result.mat, chi.sq.result.mat=chi.sq.result.mat, filenamerun=filenamerun, filenamefold=filenamefold)
}

# Goal:      Perform a single run of cross validation
# Arguments: exset = matrix of independent variables in the training set (variables = columns; samples = rows)
#            survTime = vector of survival times for all samples in the training set
#            censor = vector of censor data for all samples in the training set (0 = censored, 1 = uncensored)
#            diseaseType = string denoting the type of cancer in the training dataset (used for writing to file; default "cancer")
#	     nbest = number specifying the number of models of each size returned to bic.surv in the BMA package (default 10)
#            maxNvar = size of the active bic.surv window (default 25 - may need to be reduced with smaller training sets)
#            p = number of top-ranked genes to be included in the iterations (default 100)
#            cutPoint = threshold for separating high- from low-risk groups; enter as 50 for a 50% cutoff (default 50)
#            rp.mat = matrix for storing the random permutations
#            overall.average.results = matrix for storing high- and low-risk samples paired with censored/uncensored data across all runs
#            p.val.result.mat = matrix for holding the p-value from each fold
#	     chi.sq.result.mat = matrix for holding the chi-square statistic from each fold
#            noFolds = the number of folds in each cross validation run
#            cv = the cross validation initialization object
#            run = keeps track of what run is currently being executed (from 1:noRuns)
#            verbose = a boolean variable indicating whether or not to cat interim information to the console (default FALSE)
# Returns:   updated overall.average.results
#            updated rp.mat 
#            updated p.val.result.mat
#            updated chi.square.result.mat

crossVal.run <- function(exset, survTime, censor, diseaseType, nbest, maxNvar, p, cutPoint, rp.mat, overall.average.results, p.val.result.mat, chi.sq.result.mat, noFolds, cv, run, verbose = FALSE) { 
		
	cat (paste("******************** BEGINNING CV RUN ", run, " ********************\n", sep = ""))
	
	# Hold prediction results for this CV run
	single.run.result.mat <- mat.or.vec(2, 3)
	rownames(single.run.result.mat) <- c("Low Risk", "High Risk")
	colnames(single.run.result.mat) <- c("Censored", "Uncensored", "Percent Uncensored")

	# Generate a random permutation for current CV run
	rp <- sample(1:cv$noSamples)

	# Add permutation to the random permuation matrix 
	for (j in 1:cv$noSamples) {
		rp.mat[run, j] <- rp[j]
	}
	if (verbose == TRUE) {
		cat ("The random permutation matrix looks like this so far:\n")
		print (rp.mat[1:run ,])
	}

	# Permute the samples
	cvpass <- exset[rp ,]
	if (verbose == TRUE) {
		cat ("Here are a few samples and genes from the permuted data:\n")
		print (cvpass[1:min(5, cv$noSamples), 1:min(5, cv$noGenes)])
		cat ("Here are the row names in permuted order:\n")
		print (rownames(cvpass))
	}

	# Permute the survival time info
	survtimepass <- survTime[rp]
	if (verbose == TRUE) {
		cat ("Here are the survival times in permuted order:\n")
		print (survtimepass)
	}

	# Permute the censor info
	censorpass <- censor[rp]
	if (verbose == TRUE) {
		cat ("Here is the censor data in permuted order:\n")
		print (censorpass)
	}

	lastSample <- 0
	
	# Perform the cross validation folds for this run
	for (k in 1:noFolds) {
		cv.folds <- crossVal.fold(cvpass, survtimepass, censorpass, diseaseType, nbest, maxNvar, p, cutPoint, lastSample, noFolds, single.run.result.mat, p.val.result.mat, chi.sq.result.mat, cv, fold=k, run, verbose = verbose)
		lastSample <- cv.folds$lastSample
		single.run.result.mat <- cv.folds$single.run.result.mat
		p.val.result.mat <- cv.folds$p.val.result.mat
		chi.sq.result.mat <- cv.folds$chi.sq.result.mat
	}
	
	cat (paste("******************** END CV RUN ", run, " ********************\n", sep = ""))
	
	# print results from this run and write to file
	cat("Overall results from this run:\n")
	print (single.run.result.mat)
	write.table(paste("Results from run ", run, ":", sep = ""), cv$filenamerun, append=TRUE, quote=FALSE)
	write.table(single.run.result.mat, cv$filenamerun, append=TRUE, quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

	# Add this run's results to the overall average result matrix
	for (v in 1:2) {
		for (w in 1:2) {
			overall.average.results[v, w] <- overall.average.results[v, w] + single.run.result.mat[v, w]
		}
	}
	
	list(overall.average.results=overall.average.results, rp.mat=rp.mat, p.val.result.mat=p.val.result.mat, chi.sq.result.mat.result=chi.sq.result.mat)
}

# Goal:      Perform a single fold of a cross validation run
# Arugments: cvpass = the permuted matrix of independent variables for this fold
#            survtimepass = the permuted vector of survival times for this fold
#            censorpass = the permuted vector of censor data for this fold
#            diseaseType = string denoting the type of cancer in the training dataset (used for writing to file; default "cancer")
#	     nbest = number specifying the number of models of each size returned to bic.surv in the BMA package (default 10)
#            maxNvar = size of the active bic.surv window (default 25 - may need to be reduced with smaller training sets)
#            p = number of top-ranked genes to be included in the iterations (default 100)
#            cutPoint = threshold for separating high- from low-risk groups; enter as 50 for a 50% cutoff (default 50)
#            lastSample = number representing the last sample in the test set for this fold
#            noFolds = the number of folds in each cross validation run
#            single.run.result.mat = matrix for holding the results of a single run of cross validation
#            p.val.result.mat = matrix for holding the p-value from each fold
#	     chi.square.result.mat = matrix for holding the chi-square statistic from each fold
#            cv = the cross validation initialization object
#            fold = keeps track of what fold is currently being executed (from 1:noFolds)
#            run = keeps track of what run is currently being executed (from 1:noRuns)
#            verbose = a boolean variable indicating whether or not to cat interim information to the console (default FALSE)
# Returns:   updated lastSample
#            updated single.run.result.mat
#	     updated p.val.result.mat
#	     updated chi.square.result.mat

crossVal.fold <- function(cvpass, survtimepass, censorpass, diseaseType, nbest, maxNvar, p, cutPoint, lastSample, noFolds, single.run.result.mat, p.val.result.mat, chi.sq.result.mat, cv, fold, run, verbose = FALSE) {
	

	cat (paste("---------- BEGINNING FOLD ", fold, " ----------\n", sep = ""))

	# Pick the correct "chunk" size for this fold
	newLastSample <- floor((cv$noSamples/noFolds) * fold)
	startSample <- lastSample + 1

	# Assign training sets
	# Test set is at the very beginning of the data
	if (lastSample == 0) {
		train <- cvpass[(newLastSample + 1):cv$noSamples ,]
		trainSurv <- survtimepass[(newLastSample + 1):cv$noSamples]
		trainCens <- censorpass[(newLastSample + 1):cv$noSamples]
	}

	# Test set is at the very end of the data
	else if (newLastSample == cv$noSamples) {
		train <- cvpass[1:(startSample - 1) ,]
		trainSurv <- survtimepass[1:(startSample - 1)]
		trainCens <- censorpass[1:(startSample - 1)]
	}

	# Test set is somewhere in the middle of the data
	else {
		train <- cvpass[c(1:lastSample, (newLastSample + 1):cv$noSamples) ,]
		trainSurv <- survtimepass[c(1:lastSample, (newLastSample + 1):cv$noSamples)]
		trainCens <- censorpass[c(1:lastSample, (newLastSample + 1):cv$noSamples)]
	}

	# Assign test sets
	test <- cvpass[startSample:newLastSample ,]
	testSurv <- survtimepass[startSample:newLastSample]
	testCens <- censorpass[startSample:newLastSample]
	lastSample <- newLastSample

	# print info
	if (verbose == TRUE) {
		cat("Row names for training set:\n")
		print (rownames(train))
		cat("Row Names for test set:\n")
		print (rownames(test))
		cat("Training survival time set:\n")
		print(trainSurv)
		cat("Testing survival time set:\n")
		print (testSurv)
		cat("Training censor set:\n")
		print (trainCens)
		cat("Testing censor set:\n")
		print (testCens)
		cat("Total censor set in case you want to check the order:\n")
		print (censorpass)
	}

	# Call the iterative BMA survival algorithm
	curr.string <- paste ("_", diseaseType, "_", cv$noGenes, sep="")
	ret.bma <- iterateBMAsurv.train.predict.assess (train.dat=train, test.dat=test, surv.time.train=trainSurv, surv.time.test=testSurv, cens.vec.train=trainCens, cens.vec.test=testCens, nbest=nbest, maxNvar=maxNvar, p=p, maxIter=200000, thresProbne0=1, cutPoint=cutPoint, verbose=verbose, suff.string=curr.string)
        
        # Only calculate stats for this fold if both risk groups and both censor categories were present in test set assignments
        if(ret.bma$success==TRUE) {
		# Get results from iterative BMA survival algorithm
		results <- ret.bma$result.table
		statistics <- ret.bma$statistics

		# Calculate the p-value
		pvalue <- 1 - pchisq(statistics$chisq, 1)
		statsPerFold <- paste("Stats for fold ", fold, " in run ", run, sep="")
		statsdoc <- paste("stats_", p, "_genes_nbest_", nbest, "_maxNvar_", maxNvar, "_cutPoint_", cutPoint, ".txt", sep="")  

		# Add p-value and chi-square statistics to their corresponding matrices
		p.val.result.mat[run,fold] <- pvalue
		chi.sq.result.mat[run,fold] <- statistics$chisq

		# Write p-value and statistics to file
		write.table(statsPerFold, statsdoc, append = TRUE, quote = FALSE)
		write.table(as.matrix(statistics), statsdoc, append = TRUE, quote = FALSE)
		write.table("P-Value:", statsdoc, append = TRUE, quote = FALSE)
		write.table(pvalue, statsdoc, append = TRUE, quote = FALSE)

		# Add the current results table to the single CV run table and calculate death percentages
		cv.tab <- crossVal.tabulate(results, single.run.result.mat)
		single.run.result.mat <- cv.tab$single.run.result.mat

		# print results and write to file
		if (verbose == TRUE) {
			cat("Results from this fold:\n")
			print (results)
			cat("Accumulated results from this run so far:\n")
			print (single.run.result.mat)
		}
		write.table(paste("Results from fold ", fold, " in run ", run, ":", sep = ""), cv$filenamefold, append=TRUE, quote=FALSE)
		write.table(results, cv$filenamefold, append=TRUE, quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
		write.table("Accumulated results from this run so far:", cv$filenamefold, append=TRUE, quote=FALSE)
		write.table(single.run.result.mat, cv$filenamefold, append=TRUE, quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t") 
	}
	
	list(lastSample=lastSample, single.run.result.mat=single.run.result.mat, p.val.result.mat=p.val.result.mat, chi.sq.result.mat=chi.sq.result.mat)
}

# Goal:      Tabulate results for a single CV fold and add them to the matrix for the current CV run
# Arguments: results = table from the current fold indicating how many high- and low-risk patient samples lived and died
#                      single.run.result.mat = matrix to hold results from a single CV run
# Returns:   updated single.run.result.mat

crossVal.tabulate <- function(results, single.run.result.mat) {
	
	# Add the current results table to the single CV run table
	for (m in 1:nrow(results)) {
		for (b in 1:2) {
			if (dimnames(results)[[1]][m]==dimnames(single.run.result.mat)[[1]][b])	{
				if (ncol(results) == 2) {
					for (n in 1:2) {
						single.run.result.mat[b,n] <- single.run.result.mat[b,n] + results[m,n]
					}
				}
				else if (dimnames(results)[[2]][1] == "0") {
					single.run.result.mat[b, "Life"] <- single.run.result.mat[b, "Life"] + results[m, 1]
				}
				else if (dimnames(results)[[2]][1] == "1") {
					single.run.result.mat[b, "Death"] <- single.run.result.mat[b, "Death"] + results[m, 1]
				}
				else {
					cat("This should never happen!!!\n")
				}
			}
		}
	}

	# Calculate the death percentages
	for (q in 1:2) {
		total <- 0
		for (r in 1:2) {
			total <- total + single.run.result.mat[q, r]
		}
		single.run.result.mat[q, 3] <- 100*(single.run.result.mat[q, 2]/total)
	}
	
	list(single.run.result.mat=single.run.result.mat)
}

# Goal: Carry out final calculations for the overall average results matrix
# Arguments: overall.average.results = matrix for storing high- and low-risk samples paired with censored/uncensored data across all runs
#            noRuns = the number of cross validation runs 
# Returns:   updated overall.average.results

crossVal.final.calc <- function(overall.average.results, noRuns) {
	
	# Calculate overall average result matrix
	for (risk in 1:2) {
		for (survival in 1:2) {
			overall.average.results[risk, survival] <- overall.average.results[risk, survival]/noRuns
		}
	}
	
	# Calculate overall death percentage
	for (risk.rows in 1:2) {
		total <- 0
		for (surv in 1:2) {
			total <- total + overall.average.results[risk.rows, surv]
		}
		overall.average.results[risk.rows, 3] <- 100*(overall.average.results[risk.rows, 2]/total)
	}
	
	list(overall.average.results=overall.average.results)
}
