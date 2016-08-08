findMST2.PPI <- 
    function(object, return.MST2only=TRUE)
{
    if(!is_igraph(object)) 
        stop("'object' must be of class igraph. See package igraph for details")

    if(!is.logical(return.MST2only)) 
        stop("'return.MST2only' must be logical")

    graph.matrix <- as_adjacency_matrix(object, sparse=FALSE)
    first.mst <- mst(object)
    mst1.matrix <- as_adjacency_matrix(first.mst, sparse=FALSE)
    distmat2 <- graph.matrix - mst1.matrix
    gr2 <- graph_from_adjacency_matrix(distmat2, weighted=TRUE, mode="undirected")
    second.mst <- mst(gr2)
    mst2.matrix <- as_adjacency_matrix(second.mst, sparse=FALSE)
    MST2.matrix <- mst1.matrix + mst2.matrix
    MST2 <- graph_from_adjacency_matrix(MST2.matrix, weighted=TRUE, mode="undirected")
    if(return.MST2only) MST2 else 
        list("MST2"=MST2, "first.mst"=first.mst, "second.mst"=second.mst)
}
