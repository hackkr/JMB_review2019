
#functions:

`%ni%` = Negate(`%in%`)

#This transposes the genes expression data so that each gene is a column and each conditon is a row. 
data.t <- function(df, genelist, scaled = c("cond", "gene")) {
  num.df <- select_if(df, is.numeric)
  if (scaled == "cond") {num.df <-  scale(num.df)}
  condlist <- as.vector(colnames(num.df))
  if (nrow(num.df) == length(genelist)) {
    tmp <- as.data.frame(t(num.df), stringsAsFactors = F)
    if (scaled == "gene") {tmp <-  scale(tmp)}
    colnames(tmp) <- genelist
    rownames(tmp) <- NULL
    wide <- data.frame("cond" = condlist, tmp, stringsAsFactors = F)
    rownames(wide) <- NULL
    return(wide) 
  }
  else  {
    return("genelist wrong length")
  }
}

var.check <- function(df.long) {
  #this function will return the genes that have 0 variance. 
  #data provided should be long format that contains gene and expression columns
  genes <- df.long %>% 
    group_by(gene) %>%
    summarise(variance = var(expression, na.rm = T)) %>% 
    filter(variance == 0) %$% gene
  var.g <- df.long %>% 
    group_by(gene) %>% 
    summarise(variance = var(expression, na.rm = T)) %>%
    filter(., !is.nan(variance))
  print(paste(length(genes), "have variance of zero. The range of variances is", round(range(var.g$variance)[1], 4), "to", round(range(var.g$variance)[2], 4)))
}

var.cutoff <- function(x, percent, return.list = F) {
  #determines the variation for each gene across all conditions. Sorts by variation and returns only the genes in the top X%.
  #can also provide a absolute number of genes (greater than 1).
  genes.tot <- nrow(x)
  sorted.genes <- arrange(x, desc(var))
  if (percent > 0 && percent < 1) {
    cutoff <- ceiling(percent*genes.tot)
  } else if (percent > 1) {
    cutoff <- percent
  }
  filtered <- head(sorted.genes, cutoff)
  if (return.list == F) {
    return(min(filtered[["var"]]))
  } else {
    print(paste(nrow(filtered), "genes pass cutoff"))
    return(filtered$gene)
  }
  
}

na.cutoff <- function(x, cond, percent) {
  #calculates number of NAs allowed given a percentage of experimental conditons and return list of gene names
  ncond <- nrow(cond)
  if (percent > 0 && percent < 1) {
    cutoff <- ncond - ceiling(percent*ncond)}
  filtered <- filter(x, na.count < cutoff)
  print(paste(nrow(filtered), "genes pass cutoff.", nrow(x)-nrow(filtered), "were removed"))
  return(filtered$gene)
}

make.dend <- function(mtx.cor) {
  mtx.dist <- as.dist(1-mtx.cor)
  as.dendrogram(hclust(mtx.dist, method = "complete"))
} 

#draws all-by-all heatmap and dendrogram
heat.plus.dend <- function(corr, dend, dendcluster = 2, key = T) {
  dend <- color_branches(dend, k=dendcluster, col = plasma(dendcluster, begin = 0.1, end = 0.9))
  heatmap.2(corr, Rowv = ladderize(dend), Colv = ladderize(dend),
            dendrogram = "row", col = plasma(15, begin = 0.1, end = 0.9), 
            trace = "none", density.info = "none",
            labRow = FALSE, labCol = FALSE,
            lhei = c(1,6), lwid = c(2,6),
            key.title = NA, key.par = list(cex=1.2),
            key.xlab = "expression correlation")
}

#CALCULATES 95% CONFIDENCE INSTERVAL
conf_int95 <- function(data) {
  n <- length(data)
  error <- qt(0.975, df = n - 1) * sd(data) / sqrt(n)
  return(error)
}

#TESTS FOR FUNCTIONAL ENRICHMENT
#performs hypergeometric test on provided subset of genes/proteins relative to the genome
#code adapted from Keely Dulmage

nogtest <- function(namelist,nogfile,namecol,pvalue, cutoff = 5) {
  #namelist is a vector of protein on gene names you wnat to test for enrichment
  #nogfile is the genome-wide GETNOG output
  #p-value is significance threshold desired
  #cutoff preset prevents functional categories with less than the designated number of genes/proteins being displayed 
  
  nogs <- nogfile[nogfile[[namecol]] %in% namelist,]
  clust= table(nogs[["COG_category"]])
  resm <- matrix(0,length(clust),3) #create 0 matrix
  res <-data.frame(resm)  #make 0 df
  rownames(res) <- names(clust)
  colnames(res) <- c("probability", "expect","count")
  all=table(nogfile[["COG_category"]][nogfile[["COG_category"]] %in% nogs[["COG_category"]]])
  for (i in 1:length(clust)){   #calc expected frequencies and pval by hypergeo and append to DF
    
    res[i,1] <- phyper(clust[i], all[i], sum(all)-all[i], nrow(nogs),lower.tail=F)
    res[i,2] <- all[i]*(nrow(nogs)/sum(all))
    res[i,3] <- clust[i]
  }
  fin <- subset(res, probability <= pvalue & count >= cutoff)
  fin <- cbind("COG" = rownames(fin), fin, stringsAsFactors = F)
  row.names(fin) <- NULL
  return(fin)
}

#Splits multidomain proteins with multiple functional categories and returns a table of category sums by domain. 
#MissingasS = TRUE converts any "NA" cog.cat into "S"
all.domains <- function(x, missingasS = FALSE) {
  #x is the a freq table/DF of cog cats
  require("plyr")
  df.list <- list()
  for (i in 1:nrow(x)) {
    a <- unlist(strsplit(as.character(x$x[i]), ", "))
    b <- rep(x$freq[i], length(a))
    df <- data.frame("cogs" = a, "freq" = b)
    df.list[[i]] <- df
  }
  bind_rows(df.list) -> df
  if (missingasS == TRUE) {df$cogs[is.na(df$cogs)] <- "S" }
  df <- ddply(df, "cogs", numcolwise(sum))
  return(df)
}

#converts any "NA" cog.cat into "S"
missingasS <- function(x, colname) {
  x[[colname]][is.na(x[[colname]])] <- "S"
  sums <- ddply(x, colname, numcolwise(sum))
  return(sums)
}

#coverts frequency table to percentages
freq.as.percent <- function(x, colname) {
  x[[colname]] <- x[["freq"]]/sum(x[["freq"]])*100
  x[["freq"]] <- NULL
  return(x)
}

sig.list <-  function(x, coglist, clusterlen) {
  require("dplyr")
  df.list <- list()
  for (i in 1:length(coglist)) {
    df <- nogset(x[["query"]], x, coglist[i])
    if (i <= clusterlen) {df$cluster <- rep(1, nrow(df))} 
    else {df$cluster <- rep(2, nrow(df))}
    df.list[[i]] <- df
  }
  fin <- bind_rows(df.list)
  return(fin)
}