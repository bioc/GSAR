CVMtest <- 
    function(object, group, nperm=1000, pvalue.only=TRUE)
{
    if(!(is.matrix(object))) 
        stop("'object' must be a matrix where rows are features 
            and columns are samples")

    if(is.null(group)) 
        stop("'group' must be a vector indicating group association. 
            Possible values are 1 and 2")

    nv <- ncol(object)

    if(!is.logical(pvalue.only)) 
        stop("'pvalue.only' must be logical")

    if(length(group) != nv) 
        stop("length of 'group' must equal the number of columns in 'object'")

    if(sum(group %in% c(1,2)) < nv) 
        stop("all members in 'group' must have values 1 or 2")

    if((sum(group == 1) < 3) | (sum(group == 2) < 3)) 
        stop("there are less than 3 samples in at least one group")

    object <- object[,c(which(group == 1),which(group == 2))]
    nv1 <- sum(group == 1)
    objt <- aperm(object, c(2,1))
    Wmat <- as.matrix(dist(objt, method="euclidean", diag=TRUE, 
        upper=TRUE, p=2))
    gr <- graph_from_adjacency_matrix(Wmat, weighted=TRUE, mode="undirected")
    V(gr)[c(1:nv1)]$color <- "green"
    V(gr)[c((nv1+1):nv)]$color <- "red"
    MST <- mst(gr)
    KSranking <- HDP.ranking(MST)
    domain <- V(MST)$color
    T_perm <- array(0,c(1,nperm))

    for(itr in 1:nperm) 
    { 
        randperm <- sample(domain, replace=FALSE)
        term1 <- sum(((which(randperm[KSranking] == "green")/(nv-nv1)) - 
           (c(1:nv1)*((1/nv1)+(1/(nv-nv1)))))^2)
        term2 <- sum(((which(randperm[KSranking] == "red")/nv1) - 
           (c(1:(nv-nv1))*((1/nv1)+(1/(nv-nv1)))))^2)
        T_perm[itr] <- ((nv1*(nv-nv1)) / (nv1+(nv-nv1))^2) * (term1 + term2)
    }

    term1 <- sum(((which(domain[KSranking] == "green")/(nv-nv1)) - 
       (c(1:nv1)*((1/nv1)+(1/(nv-nv1)))))^2)
    term2 <- sum(((which(domain[KSranking] == "red")/nv1) - 
       (c(1:(nv-nv1))*((1/nv1)+(1/(nv-nv1)))))^2)
    T_obs <- ((nv1*(nv-nv1)) / (nv1+(nv-nv1))^2) * (term1 + term2)

    pvalue <- (sum(T_perm >= T_obs) + 1) / (length(T_perm) + 1)
    if(pvalue.only) return(pvalue)
    if(!pvalue.only) return(list("statistic"=T_obs,"perm.stat"=T_perm,"p.value"=pvalue))
}
