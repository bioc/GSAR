TestGeneSets <-
    function(object, group, geneSets=NULL, min.size=10, max.size=500,
        test=NULL, nperm=1000, mst.order=1, pvalue.only=TRUE)
{
    if(!(is.matrix(object))) 
        stop("'object' must be a matrix where rows are features and 
            columns are samples")

    if(is.null(group)) 
        stop("'group' must be a vector indicating group association. 
            Possible values are 1 and 2")

    if((is.null(geneSets)) | (!is.list(geneSets)))
        stop("'geneSets' must be a list object where each entry is a 
            character vector of feature (gene) identifiers")

    if((is.null(test)) | (!is.character(test)))
        stop("'test' must be one of 'GSNCAtest', 'KStest', 'MDtest', 
            'RKStest', 'RMDtest', or 'WWtest'")

    if(!is.logical(pvalue.only)) 
        stop("'pvalue.only' must be logical")

    nv <- ncol(object)

    if(length(group) != nv) 
        stop("length of 'group' must equal the number of columns in 'object'")

    if(sum(group %in% c(1,2)) < nv) 
        stop("all members in 'group' must have values 1 or 2")

    if((sum(group == 1) < 3) | (sum(group == 2) < 3)) 
        stop("there are less than 3 samples in at least one group")

    if(mst.order > 5)
    {
        warning("'mst.order' cannot be greater than 5. Larger values are 
            reduced to 5")
        mst.order <- 5
    }

    object <- object[,c(which(group == 1),which(group == 2))]
    nv1 <- sum(group == 1)
    rn <- rownames(object)
    for(k in 1:length(geneSets))
        geneSets[[k]] <- intersect(geneSets[[k]], rn)
    gs.size <- sapply(geneSets, "length")
    geneSets <- geneSets[(gs.size>=min.size) & (gs.size<=max.size)]
    res <- list()

    for(k in 1:length(geneSets))
    {
    if(test=="GSNCAtest") 
        res[[length(res)+1]] <- GSNCAtest(object=object[geneSets[[k]],],
            group=group, nperm=nperm, pvalue.only=pvalue.only)
    if(test=="WWtest") 
        res[[length(res)+1]] <- WWtest(object=object[geneSets[[k]],],
            group=group, nperm=nperm, pvalue.only=pvalue.only)
    if(test=="KStest")
        res[[length(res)+1]] <- KStest(object=object[geneSets[[k]],],
            group=group, nperm=nperm, pvalue.only=pvalue.only)
    if(test=="MDtest")
        res[[length(res)+1]] <- MDtest(object=object[geneSets[[k]],], 
            group=group, nperm=nperm, pvalue.only=pvalue.only)
    if(test=="ADtest")
        res[[length(res)+1]] <- ADtest(object=object[geneSets[[k]],],
            group=group, nperm=nperm, pvalue.only=pvalue.only)
    if(test=="CVMtest")
        res[[length(res)+1]] <- CVMtest(object=object[geneSets[[k]],],
            group=group, nperm=nperm, pvalue.only=pvalue.only)
    if(test=="RKStest")
        res[[length(res)+1]] <- RKStest(object=object[geneSets[[k]],],
            group=group, mst.order=mst.order, nperm=nperm,
                pvalue.only=pvalue.only)
    if(test=="RMDtest")
        res[[length(res)+1]] <- RMDtest(object=object[geneSets[[k]],],
            group=group, mst.order=mst.order, nperm=nperm,
                pvalue.only=pvalue.only)
    if(test=="RADtest")
        res[[length(res)+1]] <- RADtest(object=object[geneSets[[k]],],
            group=group, mst.order=mst.order, nperm=nperm,
                pvalue.only=pvalue.only)
    if(test=="RCVMtest")
        res[[length(res)+1]] <- RCVMtest(object=object[geneSets[[k]],],
            group=group, mst.order=mst.order, nperm=nperm,
                pvalue.only=pvalue.only)
    }

    if(!is.null(names(geneSets)))
        names(res) <- names(geneSets)
    return(res)
}