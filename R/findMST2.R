findMST2 <- 
    function(object, cor.method="pearson", min.sd=1e-3, 
        return.MST2only=TRUE)
{
    if(!(is.matrix(object))) 
        stop("'object' must be a matrix where rows are features and 
            columns are samples")

    if(!(cor.method %in% c("pearson", "spearman", "kendall"))) 
        stop("'cor.method' must be a character string indicating which 
            correlation coefficient to be calculated. 
                One of 'pearson' (default), 'spearman' or 'kendall'")

    if(is.null(rownames(object))) 
        rownames(object) <- as.character(c(1:nrow(object)))

    objt <- aperm(object, c(2,1))
    sdf <- apply(objt, 2, "sd")

    if(sum(sdf < min.sd) == 1) 
        stop(paste("feature ", which(sdf < min.sd), " has a standard 
            deviation smaller than 'min.sd'", sep=""))

    if(sum(sdf < min.sd) > 1) 
        stop(paste("there are ", sum(sdf < min.sd), " features with standard 
            deviation smaller than ", min.sd, sep=""))

    distmat <- 1 - abs(cor(objt, method=cor.method))
    gr <- graph_from_adjacency_matrix(distmat, weighted=TRUE, mode="undirected")
    first.mst <- mst(gr)
    mst1.matrix <- as_adjacency_matrix(first.mst, attr="weight", sparse=FALSE)
    distmat2 <- distmat - mst1.matrix
    gr2 <- graph_from_adjacency_matrix(distmat2, weighted=TRUE, mode="undirected")
    second.mst <- mst(gr2)
    mst2.matrix <- as_adjacency_matrix(second.mst, attr="weight", sparse=FALSE)
    MST2.matrix <- mst1.matrix + mst2.matrix
    MST2 <- graph_from_adjacency_matrix(MST2.matrix, weighted=TRUE, mode="undirected")
    if(return.MST2only) MST2 else 
        list("MST2"=MST2, "first.mst"=first.mst, "second.mst"=second.mst)
}
