radial.ranking <- 
    function(object)
{
    if(!is_igraph(object)) 
        stop("'object' must be of class igraph. See package igraph for details")

    mst.matrix.unweighted <- as_adjacency_matrix(object, attr=NULL, sparse=FALSE)
    mst.unweighted <- graph_from_adjacency_matrix(mst.matrix.unweighted, 
        weighted=NULL, mode="undirected")
    sp <- apply(shortest.paths(mst.unweighted), 1, max)
    spw <- apply(shortest.paths(object), 1, max)
    new.sp <- sp + spw/max(spw)
    radius <- min(new.sp)
    center <- which(new.sp == radius)
    if (length(center) > 1) center <- center[1]
    ranktree <- sort(shortest.paths(mst.unweighted)[, center] + 
        (shortest.paths(object)[, center]/max(shortest.paths(object)[, center])),
            decreasing=FALSE, index.return=TRUE)

    radial_ranking <- ranktree$ix
    return(radial_ranking)
}
