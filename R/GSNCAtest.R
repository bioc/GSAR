GSNCAtest <- 
    function(object, group, nperm=1000, cor.method="pearson",
        check.sd=TRUE, min.sd=1e-3, max.skip=10, pvalue.only=TRUE)
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

    if((sum(group==1)<3) | (sum(group==2)<3)) 
        stop("there are less than 3 samples in at least one group")

    if(!(cor.method %in% c("pearson", "spearman", "kendall"))) 
        stop("'cor.method' must be a character string indicating which 
            correlation coefficient to be calculated. 
                One of 'pearson' (default), 'spearman' or 'kendall'")

    object <- object[,c(which(group == 1), which(group == 2))]
    nv1 <- sum(group==1)
    objt <- aperm(object, c(2,1))
    group1 <- objt[1:nv1,]
    group2 <- objt[(nv1+1):nv,]

    if(check.sd == TRUE)
    {
        sd1 <- apply(group1, 2, "sd")
        sd2 <- apply(group2, 2, "sd")

        if(sum(sd1 < min.sd) == 1) 
            stop(paste("feature ", which(sd1 < min.sd), " has a standard 
                deviation smaller than 'min.sd' in group 1", sep=""))

        if(sum(sd2 < min.sd) == 1) 
            stop(paste("feature ", which(sd2 < min.sd), " has a standard 
                deviation smaller than 'min.sd' in group 2", sep=""))

        if(sum(sd1 < min.sd) > 1) 
            stop(paste("there are ", sum(sd1 < min.sd), " features with 
                standard deviation smaller than ", min.sd, " in group 1", 
                    sep=""))

        if(sum(sd2 < min.sd) > 1) 
            stop(paste("there are ", sum(sd2 < min.sd), " features with 
                standard deviation smaller than ", min.sd, " in group 2", 
                    sep=""))
    }

    cormat1 <- abs(cor(group1, method=cor.method))
    cormat2 <- abs(cor(group2, method=cor.method))
    e1 <- eigen(cormat1)
    e2 <- eigen(cormat2)
    p1 <- abs(e1$vectors[,1])
    p2 <- abs(e2$vectors[,1])
    D_obs <- sum(abs((p1*norm(matrix(p1))) - (p2*norm(matrix(p2)))))
    domain <- c(1:nv)
    D_perm <- array(0,c(1,nperm))
    skip.counter <- 0

    for(itr in 1:nperm)
    { 
        randperm <- sample(domain, replace=FALSE)
        objt <- aperm(object[,randperm], c(2,1))
        group1 <- objt[1:nv1,]
        group2 <- objt[(nv1+1):nv,]
        if(check.sd == TRUE)
        {
            sd1 <- apply(group1, 2, "sd")
            sd2 <- apply(group2, 2, "sd")
            while(((sum(sd1 < min.sd)>0) | (sum(sd2 < min.sd)>0)) & 
                (skip.counter <= max.skip))
            {
            if(skip.counter == max.skip) 
                stop("number of skipped permutations exceeded 'max.skip'")
            skip.counter <- skip.counter+1
            randperm <- sample(domain, replace=FALSE)
            objt <- aperm(object[,randperm], c(2,1))
            group1 <- objt[1:nv1,]
            group2 <- objt[(nv1+1):nv,]
            sd1 <- apply(group1, 2, "sd")
            sd2 <- apply(group2, 2, "sd")
            }
        }
        cormat1 <- abs(cor(group1, method=cor.method))
        cormat2 <- abs(cor(group2, method=cor.method))
        e1 <- eigen(cormat1)
        e2 <- eigen(cormat2)
        p1 <- abs(e1$vectors[,1])
        p2 <- abs(e2$vectors[,1])
        D_perm[itr] <- sum(abs((p1*norm(matrix(p1))) - (p2*norm(matrix(p2)))))
    }
    pvalue <- (sum(D_perm >= D_obs) + 1) / (length(D_perm) + 1)
    if(pvalue.only) return(pvalue)
    if(!pvalue.only) return(list("statistic"=D_obs,"perm.stat"=D_perm,"p.value"=pvalue))
}
