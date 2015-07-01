plotMST2.pathway <- 
    function(object, group, name=NULL, legend.size=1, 
        label.size=1, cor.method="pearson", min.sd=1e-3, return.weights=FALSE)
{
    if(!(is.matrix(object))) 
        stop("'object' must be a matrix where rows are features 
            and columns are samples")

    if(is.null(group)) 
        stop("'group' must be a vector indicating group association. 
            Possible values are 1 and 2")

    nv <- ncol(object)

    if(length(group) != nv) 
        stop("length of 'group' must equal the number of columns in 'object'")

    if(sum(group %in% c(1,2)) < nv) 
        stop("all members in 'group' must have values 1 or 2")

    if((sum(group == 1) < 3) || (sum(group == 2) < 3))
        stop("there are less than 3 samples in at least one group")

    if(!(cor.method %in% c("pearson", "spearman", "kendall"))) 
        stop("'cor.method' must be a character string indicating which 
            correlation coefficient to be calculated. One of 'pearson' 
                (default), 'spearman' or 'kendall'")

    object <- object[,c(which(group == 1), which(group == 2))]
    nv1 <- sum(group == 1)
    if(length(rownames(object)) < nrow(object)) 
        gnames <- as.character(c(1:nrow(object))) else gnames <- rownames(object)

    objt <- aperm(object, c(2,1))
    group1 <- objt[1:nv1,]
    group2 <- objt[(nv1+1):nv,]
    sd1 <- apply(group1, 2, "sd")
    sd2 <- apply(group2, 2, "sd")

    if(sum(sd1 < min.sd) == 1) 
        stop(paste("feature ", which(sd1 < min.sd), " has a standard deviation 
            smaller than 'min.sd' in group 1", sep=""))

    if(sum(sd2 < min.sd) == 1) 
        stop(paste("feature ", which(sd2 < min.sd), " has a standard deviation 
            smaller than 'min.sd' in group 2", sep=""))

    if(sum(sd1 < min.sd) > 1) 
        stop(paste("there are ", sum(sd1 < min.sd), " features with standard 
            deviation smaller than ", min.sd, " in group 1", sep=""))

    if(sum(sd2 < min.sd) > 1) 
        stop(paste("there are ", sum(sd2 < min.sd), " features with standard 
            deviation smaller than ", min.sd, " in group 2", sep=""))

    cormat1 <- abs(cor(group1, method=cor.method))
    cormat2 <- abs(cor(group2, method=cor.method))
    e1 <- eigen(cormat1)
    e2 <- eigen(cormat2)
    p1 <- matrix(abs(e1$vectors[,1]))
    p2 <- matrix(abs(e2$vectors[,1]))
    p1 <- p1 * norm(p1)
    p2 <- p2 * norm(p2)
    colnames(p1) <- "class1"
    colnames(p2) <- "class2"
    rownames(p1) <- rownames(p2) <- gnames
    major1.val <- max(p1)
    major2.val <- max(p2)
    major1.ind <- which.max(p1)
    major2.ind <- which.max(p2)
#    names(p1) <- names(p2) <- gnames
    MST2.group1 <- findMST2(object[,c(1:nv1)], cor.method, min.sd, TRUE)
    MST2.group2 <- findMST2(object[,c((nv1+1):nv)], cor.method, min.sd, TRUE)
    V(MST2.group1)$color <- V(MST2.group2)$color <- "red4"
    V(MST2.group1)$color[p1 < 1.5] <- "red"
    V(MST2.group2)$color[p2 < 1.5] <- "red"
    V(MST2.group1)$color[p1 < 1.25] <- "orange"
    V(MST2.group2)$color[p2 < 1.25] <- "orange"
    V(MST2.group1)$color[p1 < 1] <- "yellow"
    V(MST2.group2)$color[p2 < 1] <- "yellow"
    V(MST2.group1)$color[p1 < 0.75] <- "gray"
    V(MST2.group2)$color[p2 < 0.75] <- "gray"
    V(MST2.group1)$color[p1 < 0.5] <- "ghostwhite"
    V(MST2.group2)$color[p2 < 0.5] <- "ghostwhite"

    par(mfrow=c(1,2), mar=c(1,2,12,2), oma=c(1,1,4,1), cex=0.7)

    plot(MST2.group1, vertex.label.font=1, vertex.label.cex=label.size, 
        vertex.label=gnames, vertex.label.dist=1, vertex.size=8, 
            layout=layout.fruchterman.reingold, edge.width=1)

    title(paste("Group 1\n", "Hub Gene (group 1):   ", gnames[major1.ind], 
        "\n", "Weight Factor:   ", floor(1000*major1.val)/1000, "\n", 
            "Hub Gene (group 2):   ",gnames[major2.ind],"\n",
                "Weight Factor:   ",floor(1000*p1[gnames[major2.ind]])/1000, 
                    "\n","\n","\n","\n","\n", "MST2 for group 1", sep=""))

    par(xpd=NA)

    legend(x=-0.8, y=1.5, cex=legend.size, legend=c("w>1.5", "1.25<w<1.5", 
        "1<w<1.25", "0.75<w<1", "0.5<w<0.75", "w<0.5"), fill=c("red4", "red", 
            "orange", "yellow", "gray", "ghostwhite"), horiz=TRUE)

    plot(MST2.group2, vertex.label.font=1, vertex.label.cex=label.size, 
        vertex.label=gnames, vertex.label.dist=1, vertex.size=8, 
            layout=layout.fruchterman.reingold, edge.width=1)

    title(paste("Group 2\n", "Hub Gene (group 2):   ", gnames[major2.ind], 
        "\n", "Weight Factor:   ", floor(1000*major2.val)/1000, "\n", 
            "Hub Gene (group 1):   ", gnames[major1.ind],"\n",
                "Weight Factor:   ", floor(1000*p2[gnames[major1.ind]])/1000, 
                    "\n","\n","\n","\n","\n", "MST2 for group 2", sep=""))

    if(!(is.null(name))) 
        mtext(paste("Pathway: ", name, sep=""), cex=0.9, outer=TRUE, line=2)

    mtext(paste("There are ", length(gnames)," genes in this pathway", sep=""), 
        cex=0.9, outer=TRUE, line=0)

if(return.weights==TRUE) return(cbind(p1,p2))
}
