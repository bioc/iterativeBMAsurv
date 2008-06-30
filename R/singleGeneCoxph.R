# Goal:      Fit a cox PH model for each gene in the breast cancer data
#            Rank each gene using loglik (higher the better)
# Requires:  Survival library
# Arguments: trainData = matrix of independent variables in the training set (variables = columns; samples = rows)
#            survData = vector of survival times for all samples in the training set
#            censoredData = vector of censor data for all samples in the training set (0 = censored, 1 = uncensored)
# Output:    print out the loglik for each gene (sorted in descending order)
# Returns:   the sorted matrix, 1st column = sorted loglik, 2nd column = original index in the trainData

singleGeneCoxph <- function(trainData, survData, censoredData) {
	
	# Create a survival object
	mySurvObj <- Surv (survData, censoredData)  

	# Fit the coxph model for each gene and save the results
	loglikVec <- rep (0, length=ncol(trainData))
	
	for (i in 1:ncol(trainData)) {
		curr.fit <- coxph (mySurvObj ~ trainData[, i])
		loglikVec[i] <- curr.fit$loglik[2]
	}
	
	names (loglikVec) <- dimnames(trainData)[[2]]

	# Sort the log likelihood in descending order
	sortedLogLikVec <- sort (loglikVec, decreasing = TRUE, index = TRUE)

	# Results matrix to print out
	retMatrix <- matrix (0, nrow=ncol(trainData), ncol=2)
	retMatrix [, 1] <- sortedLogLikVec$x
	retMatrix [, 2] <- sortedLogLikVec$ix
	dimnames(retMatrix) <- list(names (sortedLogLikVec$x), c("sorted numbers", "original index"))
	write.table (retMatrix, file="sorted_loglik.txt", quote=FALSE, sep="\t")

	return (retMatrix)
}


# Goal: Given the sorted return matrix, print out the top G genes
# Arguments: retMatrix = the sorted return matrix
#            numGlist = list of acceptable values for top G genes
#            trainData = matrix of independent variables in the training set (variables = columns; samples = rows)
# Output:    a printed list of the top G genes

printTopGenes <- function(retMatrix, numGlist=c(10, 30, 50, 100, 500, 1000, ncol(trainData)), trainData, myPrefix="sorted_topCoxphGenes_") {

	# Get all the samples (rows), with the genes sorted 
	sorted.data <- trainData[, retMatrix[, 2]]

	# Write the sorted genes to file
	for (G in numGlist) {
		curr.data <- sorted.data[, 1:G]
		currfilename <- paste (myPrefix, G, ".txt", sep="")
		write.table (curr.data, file=currfilename, quote=FALSE, sep="\t")
	}
}

