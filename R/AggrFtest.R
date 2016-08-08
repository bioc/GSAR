AggrFtest <- 
    function(object, group, nperm=1000, pvalue.only=TRUE)
{
    if(!(is.matrix(object))) 
        stop("'object' must be a matrix where rows are features and 
            columns are samples")

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

    if((sum(group == 1) < 3) || (sum(group == 2) < 3)) 
        stop("there are less than 3 samples in at least one group")

    object <- object[,c(which(group == 1),which(group == 2))]
    ind1 <- which(group==1)
    ind2 <- which(group==2)
    p_values <- apply(object, 1, function(z) var.test(z[ind1], z[ind2])$p.value)
    p_values <- p.adjust(p_values, method="BH")
    stat_obs <- -2 * sum(log(p_values))
    stat_perm <- array(0,c(1,nperm))

    for(k in 1:nperm)
    {
        randperm <- sample(c(1:ncol(object)), replace=FALSE)
        p_values <- apply(object[,randperm], 1, function(z) var.test(z[ind1], z[ind2])$p.value)
        p_values <- p.adjust(p_values, method="BH")
        stat_perm[k] <- -2 * sum(log(p_values))
    }
    pvalue <- (sum(stat_perm >= stat_obs)+1) / (nperm+1)

    if(pvalue.only) return(pvalue)
    if(!pvalue.only) return(list("statistic"=stat_obs,"perm.stat"=stat_perm,"p.value"=pvalue))
}