library(Seurat)
library(dplyr)
library(tibble)
library(reshape2)
require(scales)
library(ggpubr)
library(tidyr)

#' 
#' @author David John
#' @param seuratObject 
#' @return filtered seurat object
FilterDeadCellsByQuantile <- function(seuratObject, lowQuantile=0.1 , highQuantile=0.95, maxMito=0.1){
  # The number of features and UMIs (nFeature_RNA and nCount_RNA) are automatically calculated for every object by Seurat.
  # For non-UMI data, nCount_RNA represents the sum of the non-normalized values within a cell
  # We calculate the percentage of mitochondrial features here and store it in object metadata as `percent.mito`.
  # We use raw count data since this represents non-transformed and non-log-normalized counts
  # The % of UMI mapping to MT-features is a common scRNA-seq QC metric.
  sizeBefore<-length(seuratObject@meta.data$orig.ident)
  cat("FilterByQuantile\n")
  #For some unknown reasons these variables need to be global for the subset function, otherwise there is an eval() unknown variable error 
  lowQuantile<<-lowQuantile
  highQuantile<<-highQuantile
  maxMito<<-maxMito
  sample<-unique(seuratObject$sample)
  Quality <- data.frame(UMI=seuratObject$nCount_RNA, nGene=seuratObject$nFeature_RNA, label = factor(seuratObject$sample), percent.mito=seuratObject$percent.mito)
  
  Quantile.low.UMI <- Quality %>% group_by(label) %>%
    summarise(UMI = list(enframe(quantile(UMI,probs = lowQuantile)))) %>%
    unnest(cols = c(UMI))
  
  Quantile.high.UMI <- Quality %>% group_by(label) %>%
    summarise(UMI = list(enframe(quantile(UMI,probs = highQuantile)))) %>%
    unnest(cols = c(UMI))
  
  Quantile.low.Gene <- Quality %>% group_by(label) %>%
    summarise(nGene = list(enframe(quantile(nGene,probs = lowQuantile)))) %>%
    unnest(cols = c(nGene))
  
  Quantile.high.Gene <- Quality %>% group_by(label) %>%
    summarise(nGene = list(enframe(quantile(nGene,probs = highQuantile)))) %>%
    unnest(cols = c(nGene))
  
  
  df<-seuratObject@meta.data
  
  gg1<- ggplot(Quality, aes(x="nUMI", y=UMI)) + geom_violin(scale = "width") + 
    theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(), legend.position = "none", axis.text.x = element_text(size=12, face = "bold"), axis.title.y = element_blank(), axis.text.y = element_text(size=10)) + 
    geom_hline(yintercept = Quantile.high.UMI$value, color="red", linetype="dashed") + geom_text(aes(0.9,Quantile.high.UMI$value, label=Quantile.high.UMI$value , vjust = -1)) + 
    geom_hline(yintercept = Quantile.low.UMI$value, color="red", linetype="dashed") + geom_text(aes(0.9,Quantile.low.UMI$value, label=Quantile.low.UMI$value , vjust = -1))
  
  gg2<- ggplot(Quality, aes(x="nFeature_RNA", y=nGene)) + geom_violin(scale = "width") + 
    theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(), legend.position = "none", axis.text.x = element_text(size=12, face = "bold"), axis.title.y = element_blank(), axis.text.y = element_text(size=10)) + 
    geom_hline(yintercept = Quantile.high.Gene$value, color="red", linetype="dashed") + geom_text(aes(0.9,Quantile.high.Gene$value, label=Quantile.high.Gene$value , vjust = -1)) +   geom_hline(yintercept = Quantile.low.Gene$value, color="red", linetype="dashed") + geom_text(aes(0.9,Quantile.low.Gene$value, label=Quantile.low.Gene$value , vjust = -1))
  
  
  gg3<- ggplot(Quality, aes(x=" % Mt Content", y=percent.mito)) + geom_violin(scale = "width") + 
    theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(), legend.position = "none", axis.text.x = element_text(size=12, face = "bold"), axis.title.y = element_blank(), axis.text.y = element_text(size=10)) + 
    geom_hline(yintercept = maxMito, color="red", linetype="dashed") + geom_text(aes(0.9,maxMito, label=maxMito , vjust = -1))
  
  gg<-ggarrange(gg1,gg2,gg3, ncol = 3)
  
  library(ggpubr)  
  
  gg<-annotate_figure(gg, fig.lab = sample, fig.lab.pos = "top", fig.lab.size = 15, fig.lab.face = "bold")
  
  seuratObject<- subset(x= seuratObject, subset = nCount_RNA < Quantile.high.UMI$value & nCount_RNA > Quantile.low.UMI$value & 
                          nFeature_RNA < Quantile.high.Gene$value & nFeature_RNA > Quantile.low.Gene$value & percent.mito < maxMito)
  
  
  
  diff<-  sizeBefore -length(seuratObject@meta.data$orig.ident)
  cat("Filtered ",diff, "from" , sizeBefore, " cells\n", "(minFeatures=",Quantile.low.Gene$value, "; maxFeatures=", Quantile.high.Gene$value, "; maxMito=" ,maxMito, ") for ", unique(seuratObject$sample), "\n" )
  rm(maxMito)
  return(list(seuratObject, gg))
}





#' Import Single cell sequencing experiments into Seurat3and perform normalisation and scale Data and do a summary of mapping stats, optional perform Doubletfinder
#' @author David John & Mariano Ruz Jurado
#' @param pathways A vector of pathways to the cellrancer count output folder (contains barcodes.tsv, genes.tsv, matrix.mtx)
#' @param ids Vector of strings that are assigned to the concordant cells
#' @return Merged seurat object
Importer <- function(pathway,id, TenX=TRUE, performNormalisation=TRUE, performScaling = FALSE,performVariableGeneDetection=TRUE, FilterCells=TRUE, FilterByAbsoluteValues=FALSE,...) {
  require(Seurat)
  if (TenX) {
    Matrix <- Read10X(pathway)
  }  else{
    Matrix <- read.table(pathway,header = TRUE,sep = ",", dec = ".", row.names = 1)
  }
  #catch optional parameters 
  optionalParameters <- list(...)
  
  seuratObject =CreateSeuratObject(counts = Matrix, project = id, min.cells = 5)
  seuratObject$sample <- id
  tmp<-unlist(strsplit(id,split = "-"))
  seuratObject$condition <- paste0(tmp[1:length(tmp)-1],collapse = "-")
  
  mito.features <- grep(pattern = "^MT-", x = rownames(x = seuratObject), value = TRUE)
  if (length(mito.features)<10) {
    mito.features <- grep(pattern = "^mt-", x = rownames(x = seuratObject), value = TRUE)
  }
  if (length(mito.features)<1) {
    warning("Error: Could not find MT genes")
    
  }
  
  percent.mito <- Matrix::colSums(x = GetAssayData(object = seuratObject, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seuratObject, slot = 'counts'))
  seuratObject$percent.mito <- percent.mito
  
  #write QC to file
  p1<-VlnPlot(object = seuratObject, features = c("nFeature_RNA"), ncol = 1, pt.size = 0) + theme(legend.position = "None", axis.title.x = element_blank(), axis.text.x = element_blank())
  p2<-VlnPlot(object = seuratObject, features = c("nCount_RNA"), ncol = 1, pt.size = 0) + theme(legend.position = "None", axis.title.x = element_blank(), axis.text.x = element_blank())
  p3<-VlnPlot(object = seuratObject, features = c("percent.mito"), ncol = 1, pt.size = 0) + theme(legend.position = "None", axis.title.x = element_blank(), axis.text.x = element_blank())
  gg_preFiltering <- ggarrange(p1,p2,p3, nrow = 1)
  gg_preFiltering <- annotate_figure(gg_preFiltering, top = text_grob(id,face="bold",color = "darkred",size=18,hjust = 0.2))
  ggsave(filename = paste0(pathway,"QC_preFiltered.svg"),device = "svg", width = 10,height = 10)
  
  if (FilterCells==TRUE) {
    message("start Filtering")
    if (FilterByAbsoluteValues==TRUE) {
      if (is.null(optionalParameters$minFeatures)) {
        stop("Please define 'minFeatures' while filtering for absolute values (FilterByAbsoluteValues==TRUE)")
      }
      if (is.null(optionalParameters$maxFeatures)) {
        stop("Please define 'maxFeatures' while filtering for absolute values (FilterByAbsoluteValues==TRUE)")
      }
      if (is.null(optionalParameters$minCounts)) {
        stop("Please define 'minCounts' while filtering for absolute values (FilterByAbsoluteValues==TRUE)")
      }
      if (is.null(optionalParameters$maxCounts)) {
        stop("Please define 'maxCounts' while filtering for absolute values (FilterByAbsoluteValues==TRUE)")
      }
      if (is.null(optionalParameters$maxMito)) {
        stop("Please define 'maxMito' while filtering for absolute values (FilterByAbsoluteValues==TRUE)")
      }
      message("Running FilterDeadCells")
      seuratObject<-FilterDeadCells(seuratObject = seuratObject,
                                    minFeatures = optionalParameters$minFeatures,
                                    maxFeatures = optionalParameters$maxFeatures,
                                    minCounts = optionalParameters$minCounts,
                                    maxCounts = optionalParameters$maxCounts,
                                    maxMito = optionalParameters$maxMito)
    }
    else {
      tmp<-FilterDeadCellsByQuantile(seuratObject = seuratObject, lowQuantile = 0.1, highQuantile = 0.95)
      seuratObject<-tmp[[1]]
      svg(paste0(pathway,"QC_QuantileFiltering.svg"))
      print(tmp[[2]])
      dev.off()
      gg_preFiltering<-tmp[[2]]
      
    }
    
  }
  if (performNormalisation==TRUE) {
    seuratObject<-NormalizeData(object = seuratObject,verbose = FALSE)
  }
  if(performVariableGeneDetection==TRUE){
    seuratObject<-FindVariableFeatures(object = seuratObject, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  }
  if (performScaling==TRUE) {
    seuratObject<-ScaleData(object = seuratObject)
  }
  message("Imported ", length(seuratObject@meta.data$orig.ident), " cells from ", pathway, "with ID ", id, "\n")
  
  
  return(list(seuratObject, gg_preFiltering))
}

#' @author Mariano Ruz Jurado
#' @param pathway A vector of pathways to the cellrancer count output folder (contains barcodes.tsv, genes.tsv, matrix.mtx)
#' @param id Vector of strings that are assigned to the concordant cells
#' @return a list with summary about mapping results 
SummarizeMapping <- function(pathway,id){
  ## mapping stats, i hope everyone uses cellranger and starsolo directories for these analysis, else no summary
  if (file.exists(paste0(unlist(strsplit(pathway,split = "outs"))[1],"outs/metrics_summary.csv"))) {
    metrics_summ.path <- paste0(unlist(strsplit(pathway,split = "outs"))[1],"outs/metrics_summary.csv")
    #define your own numeric class
    setClass('myNum')
    #define conversion
    setAs("character", "myNum",
          function(from) as.numeric(gsub(",","", gsub("%","",from))))
    #read data with custom colClasses
    metrics_summ <- read.csv(metrics_summ.path,
                             header = T,   
                             stringsAsFactors=FALSE,
                             colClasses=c("myNum"))
    
    typeof(metrics_summ$Fraction.Reads.in.Cells)
    
    metrics_col <- as.data.frame(colnames(metrics_summ))
    rownames(metrics_col) <- metrics_col[,1]
    metrics_col[,1] <- as.character(as.vector(metrics_summ[1,]))
    metrics_summ <- metrics_col
    
    # warnings CELLRANGER
    if (metrics_summ[grep(pattern = "Confidently.to.Genome",rownames(metrics_summ)),] < 70) {
      warning(paste0(id,": Reads mapped confidently to genome only ", metrics_summ[grep(pattern = "Confidently.to.Genome",rownames(metrics_summ)),]))
    }
    if (metrics_summ[grep(pattern = "Confidently.to.Transcriptome",rownames(metrics_summ)),] < 30) {
      warning(paste0(id,": Reads mapped confidently to transcriptome only ", metrics_summ[grep(pattern = "Confidently.to.Transcriptome",rownames(metrics_summ)),]))
    }
    if (paste(unlist(strsplit(metrics_summ[grep(pattern = "Number.of.Cells",rownames(metrics_summ)),],split=",")),collapse = "") < 1000) {
      warning(paste0(id,": Estimated Number of Cells only ", metrics_summ[grep(pattern = "Number.of.Cells",rownames(metrics_summ)),]), " ,maybe the 0s were cut because of CR way of displaying numbers,  if unsure check CR web_summary")
    }
    if (as.numeric(paste(unlist(strsplit(metrics_summ[grep(pattern = "Median.Genes.per.Cell",rownames(metrics_summ)),],split=",")),collapse = "")) < 300) {
      warning(paste0(id,": Median Genes per Cell only ", metrics_summ[grep(pattern = "Median.Genes.per.Cell",rownames(metrics_summ)),])," ,maybe the 0s were cut because of CR way of displaying numbers, if unsure check CR web_summary")
    }
  } else {
    metrics_summ.path <- paste0(unlist(strsplit(pathway,split = "Gene"))[1],"Gene/Summary.csv")
    metrics_summ <- read.delim2(metrics_summ.path, header = F, sep = ",")
    rownames(metrics_summ) <- metrics_summ[,1]
    metrics_summ[,1] <- NULL
    
    # warnings STAR
    if (metrics_summ[7,] < 0.70) { # mapped to genome, no grep since same name as other row 
      warning(paste0(id,": Reads mapped confidently to genome only ",metrics_summ[7,]))
    }
    if (metrics_summ[grep(pattern = "Transcriptome: Unique Genes",rownames(metrics_summ)),] < 0.30) {
      warning(paste0(id,": Reads mapped confidently to transcriptome only ", metrics_summ[grep(pattern = "Transcriptome: Unique Genes",rownames(metrics_summ)),]))
    }
    if (metrics_summ[grep(pattern = "Number of Cells",rownames(metrics_summ)),] < 1000) {
      warning(paste0(id,": Estimated Number of Cells only ", metrics_summ[grep(pattern = "Number of Cells",rownames(metrics_summ)),]))
    }
    if (metrics_summ[grep(pattern = "Median Genes per Cell",rownames(metrics_summ)),] < 300) {
      warning(paste0(id,": Median Genes per Cell only ", metrics_summ[grep(pattern = "Median Genes per Cell",rownames(metrics_summ)),]))
    }
  } 
  return(metrics_summ)
}

#' 
#' @author David John & Mariano Ruz Jurado
#' @param seuratObject 
#' @return filtered seurat object
FilterDeadCells <- function(seuratObject, minFeatures=300, maxFeatures=6000,minCounts=500,maxCounts=15000, maxMito=0.05){
  # The number of features and UMIs (nFeature_RNA and nCount_RNA) are automatically calculated for every object by Seurat.
  # For non-UMI data, nCount_RNA represents the sum of the non-normalized values within a cell
  # We calculate the percentage of mitochondrial features here and store it in object metadata as `percent.mito`.
  # We use raw count data since this represents non-transformed and non-log-normalized counts
  # The % of UMI mapping to MT-features is a common scRNA-seq QC metric.
  sizeBefore<-length(seuratObject@meta.data$orig.ident)
  
  #For some unknown reasons these variables need to be global for the subset function, otherwise there is an eval() unknown variable error 
  minFeatures<<-minFeatures
  maxFeatures<<-maxFeatures
  maxMito<<-maxMito
  seuratObject <- subset(x = seuratObject, subset = nFeature_RNA > minFeatures & nFeature_RNA < maxFeatures & nCount_RNA > minCounts & nCount_RNA < maxCounts & percent.mito < maxMito)
  
  diff<-  sizeBefore -length(seuratObject@meta.data$orig.ident)
  percent <- round(diff/sizeBefore*100,digits = 2)
  cat("Filtered ",diff, "from" , sizeBefore, " cells [",percent,"%]\n", "(minFeatures=",minFeatures, "; maxFeatures=", maxFeatures, "; minCounts=" ,minCounts,  "; maxCounts=" ,maxCounts , "; maxMito=" ,maxMito, ") for ", unique(seuratObject$sample), "\n" )
  rm(minFeatures,maxFeatures,maxMito)
  return(seuratObject)
}


## prepare data for cluster t.test from the deg list and do a cluster t-test
do_cluster_t_test <- function(seurat_subset, degs, group="condition", cluster="seurat_clusters"){
  gene_names<- names(table(rownames(degs)))
  #print(head(gene_names))
  p_values <- vector("list",length(gene_names))
  names(p_values) <- gene_names
  #gene_names <- row.names(cluster_subset)
  #if (celltype=="Adipocytes"){
  #  seurat_subset <- seurat_subset[,!seurat_subset$orig.ident=="D7"]
  #}
  group <- seurat_subset[[group]][,1]
  cluster <- seurat_subset[[cluster]][,1]
  for (gene in gene_names){
    y <- c(t(as.matrix(seurat_subset@assays$RNA[gene,])))
    test_info <- my.t.test.cluster(y, cluster, group)
    p_values[[gene]] <- test_info[nrow(test_info)]
  }
  return(p_values)
}

## added line 54-56 so that each group is tested if 
## only one observation is present and throw an error
my.t.test.cluster <- function (y, cluster, group, conf.int = 0.95) 
{
  group <- as.factor(group)
  cluster <- as.factor(cluster)
  s <- !(is.na(y) | is.na(cluster) | is.na(group))
  y <- y[s]
  cluster <- cluster[s]
  group <- group[s]
  n <- length(y)
  if (n < 2) 
    stop("n<2")
  gr <- levels(group)
  if (length(gr) != 2) 
    stop("must have exactly two treatment groups")
  n <- table(group)
  nc <- tapply(cluster, group, function(x) length(unique(x)))
  bar <- tapply(y, group, mean)
  u <- unclass(group)
  y1 <- y[u == 1]
  y2 <- y[u == 2]
  c1 <- factor(cluster[u == 1])
  c2 <- factor(cluster[u == 2])
  b1 <- tapply(y1, c1, mean)
  b2 <- tapply(y2, c2, mean)
  m1 <- table(c1)
  m2 <- table(c2)
  if (any(names(m1) != names(b1)))
    stop("logic error 1")
  if (any(names(m2) != names(b2)))
    stop("logic error 2")
  if (any(m2 < 2))
    stop(paste("The following clusters contain only one observation:",
               paste(names(m2[m2 < 2]), collapse = " ")))
  if (any(m1 < 2))
    stop(paste("The following clusters contain only one observation:",
               paste(names(m1[m1 < 2]), collapse = " ")))
  M1 <- mean(y1)
  M2 <- mean(y2)
  ssc1 <- sum(m1 * ((b1 - M1)^2))
  ssc2 <- sum(m2 * ((b2 - M2)^2))
  if (nc[1] != length(m1))
    stop("logic error 3")
  if (nc[2] != length(m2))
    stop("logic error 4")
  df.msc <- sum(nc) - 2
  msc <- (ssc1 + ssc2)/df.msc
  v1 <- tapply(y1, c1, var)
  v2 <- tapply(y2, c2, var)
  ssw1 <- sum((m1 - 1) * v1)
  ssw2 <- sum((m2 - 1) * v2)
  df.mse <- sum(n) - sum(nc)
  mse <- (ssw1 + ssw2)/df.mse
  na <- (sum(n) - (sum(m1^2)/n[1] + sum(m2^2)/n[2]))/(sum(nc) - 
                                                        1)
  rho <- (msc - mse)/(msc + (na - 1) * mse)
  r <- max(rho, 0)
  C1 <- sum(m1 * (1 + (m1 - 1) * r))/n[1]
  C2 <- sum(m2 * (1 + (m2 - 1) * r))/n[2]
  v <- mse * (C1/n[1] + C2/n[2])
  v.unadj <- mse * (1/n[1] + 1/n[2])
  de <- v/v.unadj
  dif <- diff(bar)
  se <- sqrt(v)
  zcrit <- qnorm((1 + conf.int)/2)
  cl <- c(dif - zcrit * se, dif + zcrit * se)
  z <- dif/se
  P <- 2 * pnorm(-abs(z))
  stats <- matrix(NA, nrow = 20, ncol = 2, dimnames = list(c("N", 
                                                             "Clusters", "Mean", "SS among clusters within groups", 
                                                             "SS within clusters within groups", "MS among clusters within groups", 
                                                             "d.f.", "MS within clusters within groups", "d.f.", "Na", 
                                                             "Intracluster correlation", "Variance Correction Factor", 
                                                             "Variance of effect", "Variance without cluster adjustment", 
                                                             "Design Effect", "Effect (Difference in Means)", "S.E. of Effect", 
                                                             paste(format(conf.int), "Confidence limits"), "Z Statistic", 
                                                             "2-sided P Value"), gr))
  stats[1, ] <- n
  stats[2, ] <- nc
  stats[3, ] <- bar
  stats[4, ] <- c(ssc1, ssc2)
  stats[5, ] <- c(ssw1, ssw2)
  stats[6, 1] <- msc
  stats[7, 1] <- df.msc
  stats[8, 1] <- mse
  stats[9, 1] <- df.mse
  stats[10, 1] <- na
  stats[11, 1] <- rho
  stats[12, ] <- c(C1, C2)
  stats[13, 1] <- v
  stats[14, 1] <- v.unadj
  stats[15, 1] <- de
  stats[16, 1] <- dif
  stats[17, 1] <- se
  stats[18, ] <- cl
  stats[19, 1] <- z
  stats[20, 1] <- P
  attr(stats, "class") <- "t.test.cluster"
  stats
}

#perform SEM Graphs
#' @author Mariano Ruz Jurado
#' @param SeuratObject # combined object
#' @param Features # vector containing featurenames
#' @param ListTest # List for which conditions t-test will be performed, if NULL always against provided CTRL 
#' @param returnValues # return df.melt.sum data frame containing means and SEM for the set group
#' @param ctrl.condition # set your ctrl condition, relevant if running with empty comparison List
#' @param group.by # select the seurat object slot where your conditions can be found, default conditon
DO.Mean.SEM.Graphs <- function(SeuratObject, Features, ListTest=NULL, returnValues=FALSE, ctrl.condition=NULL, group.by = "condition"){ 
  require(ggplot2)
  require(ggpubr)
  #require(tidyverse)
  require(reshape2)
  #SEM function definition
  SEM <- function(x) sqrt(var(x)/length(x))
  #create data frame with conditions from provided SeuratObject, aswell as original identifier of samples
  df<-data.frame(condition=setNames(as.character(SeuratObject[[group.by]][,group.by]), rownames(SeuratObject[[group.by]]))
                 ,orig.ident = SeuratObject$orig.ident)
  #get expression values for genes from individual cells, add to df
  for(i in Features){
    df[,i] <- SeuratObject@assays$RNA@data[i,]
  }
  
  #melt results 
  df.melt <- melt(df)
  #group results and summarize, also add/use SEM 
  df.melt.sum <- df.melt %>% 
    group_by(condition, variable) %>% 
    summarise(Mean = mean(value), SEM = SEM(value))
  #second dataframe containing mean values for individual samples
  df.melt.orig <- df.melt %>% 
    group_by(condition, variable, orig.ident) %>% 
    summarise(Mean = mean(value))
  
  #create comparison list for t.test, always against control, so please check your sample ordering
  # ,alternative add your own list as argument
  if (is.null(ListTest)) {
    #if ListTest is empty, so grep the ctrl conditions out of the list 
    # and define ListTest comparing every other condition with that ctrl condition
    cat("ListTest empty, comparing every sample with provided control")
    conditions <- unique(SeuratObject[[group.by]][,group.by])
    #set automatically ctrl condition if not provided
    if (is.null(ctrl.condition)) { 
      ctrl.condition <- conditions[grep(pattern = paste(c("CTRL","Ctrl","WT","Wt","wt"),collapse ="|")
                                        ,conditions)[1]]
    }
    
    df.melt.sum$condition <- factor(df.melt.sum$condition
                                    ,levels = c(ctrl.condition,levels(factor(df.melt.sum$condition))[!(levels(factor(df.melt.sum$condition)) %in% ctrl.condition)]))
    #create ListTest
    ListTest <- list()
    for (i in 1:length(conditions)) {
      cndtn <- conditions[i] 
      if(cndtn!=ctrl.condition)
      {
        ListTest[[i]] <- c(ctrl.condition,cndtn)
      }
    }
  }
  #delete Null values, created by count index
  ListTest <- ListTest[!sapply(ListTest, is.null)]
  #create barplot with significance
  p<-ggplot(df.melt.sum, aes(x = condition, y = Mean, fill = condition))+
    geom_col(color = "black")+
    geom_errorbar(aes(ymin = Mean-SEM, ymax = Mean+SEM), width = .2)+
    geom_jitter(data = df.melt.orig, aes(x=condition,y=Mean), size = 1, shape=1)+
    #ordering, control always first
    scale_x_discrete(limits=c(ctrl.condition,levels(factor(df.melt.sum$condition))[!(levels(factor(df.melt.sum$condition)) %in% ctrl.condition)]))+
    #t-test, always against control, using means from orig sample identifier
    stat_compare_means(data=df.melt.orig, comparisons = ListTest, method = "t.test", size=3)+
    facet_wrap(~variable, ncol = 9, scales = "free") +
    scale_fill_manual(values = rep(c("#FFFFFF","#BDD7E7" ,"#6BAED6", "#3182BD", "#08519C"),5) #20 colours set for more change number
                      , name = "Condition")+
    labs( title = "", y = "Mean UMI") +
    theme_classic() +
    theme(axis.text.x = element_text(face = "bold", color = "black",angle = 45,hjust = 1, size = 18),
          axis.text.y = element_text(color = "black", size = 18),
          axis.title.x = element_blank(),
          axis.title = element_text(size = 18, color = "black"),
          plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
          axis.line = element_line(color = "black", size = .6),
          strip.text.x = element_text(size = 25, color = "black"),
          legend.position = "bottom")
  print(p)
  if (returnValues==TRUE) {
    return(df.melt.sum)
  }
}

# GO analysis over-representation analysis for the subontologies Biological Process, Molecular function and Cellular Component, works for mice and human
#' @author Mariano Ruz Jurado
#' @param Differential.expression.feature.list A list containing sub-lists with names for species and their expression features
#' @param SeuratObject Single Cell Object created by Seurat
#' @param pvalueCutoff p-value Cutoff
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param qvalueCutoff qvalue cutoff on enrichment tests to report as significant. Tests must pass i) pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff on adjusted pvalues and iii) qvalueCutoff on qvalues to be reported.
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize maximal size of genes annotated for testing
#' @return list with all GO results
GO.Analysis <- function(DEG.list,SeuratObject,pvalueCutoff,pAdjustMethod,qvalueCutoff,minGSSize,maxGSSize){
  require(org.Hs.eg.db)
  require(org.Mm.eg.db)
  require(clusterProfiler)
  GO.result <- list()
  SeuratObject.combined <- SeuratObject
  for (i in 1:length(names(DEG.list))) {
    
    i=names(DEG.list[i])
    species=names(DEG.list[i])
    print(species)
    
    if ("Human" %in% unlist(strsplit(species, split = "_")) && !("integrated" %in% unlist(strsplit(species, split = "_")))) {
      OrgDb.name <- "org.Hs.eg.db"
      OrgDb = org.Hs.eg.db
      SeuratObject.universe <- SeuratObject.combined
    }
    
    if ("Mice" %in% unlist(strsplit(species, split = "_")) && !("integrated" %in% unlist(strsplit(species, split = "_")))){
      OrgDb.name <- "org.Hs.eg.db"
      OrgDb = org.Hs.eg.db
      SeuratObject.universe <- SeuratObject.combined
    }
    
    # if (c("Human","integrated") %in% unlist(strsplit(species, split = "_"))) {
    #   OrgDb.name <- "org.Hs.eg.db"
    #   OrgDb = org.Hs.eg.db
    #   SeuratObject.universe <- SeuratObject.combined
    # }
    
    if (!("Human" %in% unlist(strsplit(species, split = "_"))) && !("Mice" %in% unlist(strsplit(species, split = "_")))) {
      stop("species not supported")
      
    }
    
    for (j in names(DEG.list[[species]])) {
      j=names(DEG.list[[species]][j])
      print(j)
      marker.contrast <- j
      DEGFolder <- paste0(outputFolder,"/",species,".",j)
      
      if(!dir.exists(DEGFolder)) {
        dir.create(DEGFolder)
      } else{"Directory exists!"}
      
      # define universe with ENTREZIDs as data frames
      egallgenes <- as.vector(rownames(SeuratObject.universe)) 
      egallgenes = bitr(egallgenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb.name)
      head(egallgenes)
      
      # overmit differentially regulated genes with ENTREZIDs as data frame
      egdiffgenes <- as.vector(rownames(DEG.list[[species]][[j]]))
      egdiffgenes = bitr(egdiffgenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb.name)
      head(egdiffgenes)
      write.csv(egdiffgenes, file = paste0(DEGFolder, "/",marker.contrast,"_","DEGs.csv"))
      # enrichGO on all subontologies for differentially regulated genes
      for (m in c("BP", "MF", "CC")) {
        print(m)
        
        #define geneList (ENTREZ IDs and LOG2FC)
        DFgenes <- DEG.list[[species]][[j]]
        
        
        egallgenes$stat <- DFgenes[["avg_log2FC"]][match(egallgenes$SYMBOL,rownames(DFgenes))] #get log2FC for DEGs in egdiffgenes
        mygeneList <- egallgenes$stat
        names(mygeneList) <- egallgenes$ENTREZID
        mygeneList<-sort(mygeneList,decreasing = T)
        mygeneList <- mygeneList[!is.na(mygeneList)]
        head(mygeneList)
        length(mygeneList)
        
        # # top 500 from abs values
        # if(length(mygeneList) > 500){
        #   abs.list <- abs(mygeneList)
        #   abs.list <- sort(abs.list, decreasing = T)
        #   gene <- mygeneList[names(abs.list[1:500])]
        # } else{
        #   gene <- mygeneList
        # }
        
        ego <- enrichGO(gene = names(mygeneList),
                        OrgDb = OrgDb,
                        keyType = "ENTREZID",
                        ont = m,
                        universe = egallgenes[,2],
                        pvalueCutoff = pvalueCutoff, 
                        pAdjustMethod = pAdjustMethod, 
                        qvalueCutoff = qvalueCutoff, 
                        minGSSize = minGSSize,
                        maxGSSize = maxGSSize,
                        readable = TRUE,
                        pool = FALSE)
        
        
        ego <- setReadable(ego, OrgDb = OrgDb, keyType = "ENTREZID")
        openxlsx::write.xlsx(ego@result,file = paste0(DEGFolder,"/GO_",m,".xlsx"))
        GO.result[[species]][[paste0(marker.contrast,"_",m,"_result")]] <- ego@result
        
        
        word.split <- strsplit(ego@result$Description," ")
        
        # head(ego@result$Description)
        head(word.split)
        ### linebreaks for plot, so the name doesnt become too long for the individual terms
        for (b in 1:length(word.split)) 
        {
          
          for (d in 1:length(word.split[[b]])) { ### checking for hyphen in words making them too long, erase them
            if (length(unlist(strsplit(word.split[[b]], " "))) < length(unlist(strsplit(word.split[[b]], "-"))))
            {
              nohyphen <- unlist(strsplit(word.split[[b]], "-"))
              word.split[[b]] <- nohyphen
            }
          }
          
          for (c in 1:length(word.split[[b]])) {  ### checking word length
            if (length(unlist(strsplit(word.split[[b]][c], split =""))) > 10)
            {
              word.split2 <- unlist(strsplit(word.split[[b]][c], split =""))
              word.split2 <- word.split2[1:10]
              word.split2 <- paste(word.split2, collapse = "")
              word.split[[b]][c] <- word.split2
            }
          }
          
          
          
          if (length(word.split[[b]]) > 4) ### linebreaks
          {
            linebreaks <- word.split[[b]][1:4]
            linebreaks <- paste(linebreaks, collapse = " ")
            linebreaks2 <- word.split[[b]][5:length(word.split[[b]])] 
            if (length(linebreaks2) > 4) 
            { 
              linebreaks2 <- linebreaks2[1:4]
              
              
            }
            linebreaks2 <- paste(linebreaks2, collapse = " ")
            word.split[[b]] <- paste(linebreaks,linebreaks2, sep = "\n" )
            
          } else (paste(word.split[[b]], collapse = " ") -> word.split[[b]])
          
        }
        ego@result$Description <- word.split
        
        GO_File <- ego
        GO_File.Data.Frame <- as.data.frame(GO_File)
        
        # Calculating Gene Ratio so it is numeric and usable for a bubbleplot
        if (length(GO_File.Data.Frame$GeneRatio != 0)) {
          for (i in 1:length(GO_File.Data.Frame$GeneRatio)){
            
            Gene.Ratio <- GO_File.Data.Frame$GeneRatio[i] 
            Gene.Ratio.calculated <- as.numeric(unlist(strsplit(Gene.Ratio, split = "/" ))[1]) / as.numeric(unlist(strsplit(Gene.Ratio, split = "/" ))[2])
            GO_File.Data.Frame$GeneRatio[i] <- Gene.Ratio.calculated
          }
          GO_File.Data.Frame$GeneRatio <- as.numeric(GO_File.Data.Frame$GeneRatio)
          
          
          
          
          
          # Bubble plot for first 10 results or for all results if less 20
          if (length(GO_File.Data.Frame$Description) > 20) {
            Set_Size <- GO_File.Data.Frame[1:20,]
          } else (Set_Size <- GO_File.Data.Frame[1:length(GO_File.Data.Frame$Description),])
          
          #Bubble plot
          pdf(paste0(DEGFolder, "/GO_Bubbleplot_",m,".pdf"))
          ggplot(Set_Size, aes(x = GeneRatio , y = fct_reorder(unlist(Description), GeneRatio ))) +
            geom_point(aes(size = GeneRatio * 455, color = p.adjust)) + 
            labs(size = "Counts") +
            theme_bw(base_size = 14) +
            theme(axis.text.y=element_text(size = 16), axis.text.x = element_text(size = 14))+
            scale_colour_gradient(limits=c(min(GO_File.Data.Frame$p.adjust), max(GO_File.Data.Frame$p.adjust)), low="blue") +
            ylab(NULL) +
            expand_limits(x = c(min(Set_Size$GeneRatio)-(min(Set_Size$GeneRatio)/100*15), max(Set_Size$GeneRatio)+(max(Set_Size$GeneRatio)/100*15))) +
            ggtitle(paste0(m))-> x
          x + guides(size = guide_legend(order = 1)) -> x
          print(x)
          dev.off()
          print(x)
          ggsave(filename = paste0(DEGFolder, "/GO_Bubbleplot_",m,".svg"), width = 10, height = 10)
          # barplot
          GO_File <- ego
          GO_File.Data.Frame <- as.data.frame(GO_File)
          
          if (length(rownames(GO_File.Data.Frame)) < 1) {  ### if data frame is empty
            warning("No significant terms detected")
          } else {
            for (i in 1:length(GO_File.Data.Frame$GeneRatio)){
              GO_File.Data.Frame$GeneRatio[i] <- as.numeric(unlist(strsplit(GO_File.Data.Frame$GeneRatio[i], split = "/" ))[1]) /
                as.numeric(unlist(strsplit(GO_File.Data.Frame$GeneRatio[i], split = "/" ))[2])
            }
            GO_File.Data.Frame$GeneRatio <- as.numeric(GO_File.Data.Frame$GeneRatio)
            
            GO_File.Data.Frame2 <- GO_File.Data.Frame
            GO_File.Data.Frame2 <- GO_File.Data.Frame2[order(GO_File.Data.Frame2$p.adjust, decreasing = F),]
            GO_File.Data.Frame <-GO_File.Data.Frame2
            
            GO_File.Data.Frame$Description <- factor(GO_File.Data.Frame$Description, levels = unique(as.character(GO_File.Data.Frame$Description)))
            GO_File.Data.Frame <- transform(GO_File.Data.Frame, Description = reorder(Description, -p.adjust))
            if (length(rownames(GO_File.Data.Frame)) < 20) {k =length(rownames(GO_File.Data.Frame))} else {k=20}
            pdf(paste0(DEGFolder, "/GO_Barplot_",m,".pdf"))
            ggplot(GO_File.Data.Frame[1:k,], aes(reorder(GO_File.Data.Frame, p.adjust), x = Description , y = -log10(p.adjust))) +
              coord_flip() +
              geom_bar(aes( fill = Count), stat="identity") +
              labs(size = "Counts") +
              theme_bw(base_size = 14) +
              theme(axis.text.y=element_text(size = 16)) +
              scale_colour_gradient(limits=c(min(GO_File.Data.Frame$Count), max(GO_File.Data.Frame$Count)), low="blue") +
              ylab("-log10(p.adjusted)") -> x
            x + guides(size = guide_legend(order = 1)) -> x
            print(x)
            dev.off()
            print(x)
            ggsave(filename = paste0(DEGFolder, "/GO_Barplot_",m,".svg"), width = 10, height = 10)
            # gene concept plot
            if (length(rownames(GO_File.Data.Frame)) >= 5) {
              pdf(paste0(DEGFolder, "/GO_GeneConcept_",m,".pdf"))          
              cnet.plot <- cnetplot(ego, foldChange=mygeneList, node_label="all",cex_label_category = 1 , cex_label_gene = 0.6)
              print(cnet.plot)
              dev.off()
              print(cnet.plot)
              ggsave(filename = paste0(DEGFolder, "/GO_GeneConcept_",m,".svg"), width = 16, height = 10)
              # circular gene concept plot
              pdf(paste0(DEGFolder, "/GO_GeneConcept_circular_",m,".pdf")) 
              cnet.plot.circ <- cnetplot(ego, foldChange = mygeneList, circular = TRUE, colorEdge = TRUE,cex_label_category = 1 , cex_label_gene = 0.6, showCategory = 4)
              print(cnet.plot.circ)
              dev.off()
              print(cnet.plot.circ)
              ggsave(filename = paste0(DEGFolder, "/GO_GeneConcept_circular_",m,".svg"), width = 16, height = 10)
            }
          }
        }
      }  
    }
  }
  return(GO.result)  
}

# Go Gene set enrichment analysis using clusterProfiler
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler)
#' @author Mariano Ruz Jurado
#' @param Gene.list A list containing sub-lists with names for species and contrast, containing genes with avg_log2FC (e.g. Gene.list$HFrEF_Human_Cardiomyocytes$Markers$avg_log2FC)
#' @param SeuratObject Single Cell Object created by Seurat
#' @param pvalueCutoff p-value Cutoff
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param qvalueCutoff qvalue cutoff on enrichment tests to report as significant. Tests must pass i) pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff on adjusted pvalues and iii) qvalueCutoff on qvalues to be reported.
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize maximal size of genes annotated for testing
#' @return list with all GSEA results
GSEA.Analysis <- function(Gene.list,SeuratObject,pvalueCutoff,pAdjustMethod,qvalueCutoff,minGSSize,maxGSSize){
  require(org.Hs.eg.db)
  require(org.Mm.eg.db)
  require(clusterProfiler)
  require(ggnewscale)
  GSEA.result <- list()
  for (i in 1:length(names(Gene.list))) {
    
    i=names(Gene.list[i]) # get naming from list for automated plot naming
    species=names(Gene.list[i]) 
    print(species)
    
    if ("Human" %in% unlist(strsplit(species, split = "_"))) {
      OrgDb.name <- "org.Hs.eg.db"
      OrgDb = org.Hs.eg.db
      SeuratObject.universe <- SeuratObject
    }
    
    if ("Mice" %in% unlist(strsplit(species, split = "_"))) { #Integration changed nomenclature to human gene nomenclature, therefore "org.Hs.eg.db"
      OrgDb.name <- "org.Hs.eg.db"
      OrgDb = org.Hs.eg.db
      SeuratObject.universe <- SeuratObject
    }
    
    
    if (!("Human" %in% unlist(strsplit(species, split = "_"))) && !("Mice" %in% unlist(strsplit(species, split = "_")))) {
      stop("species not supported")
    }
    
    # define universe with ENTREZIDs as data frames
    egallgenes <- as.vector(rownames(SeuratObject.universe)) 
    egallgenes = bitr(egallgenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb.name)
    head(egallgenes)    
    
    
    
    for (j in names(Gene.list[[species]])) {
      j=names(Gene.list[[species]][j])
      marker.contrast <- j
      DEGFolder <- paste0(outputFolder,"/",species,".",j)
      
      if(!dir.exists(DEGFolder)) {
        dir.create(DEGFolder)
      } else{"Directory exists!"}
      
      DFgenes <- Gene.list[[species]][[j]]
      
      
      egallgenes$stat <- DFgenes[["avg_log2FC"]][match(egallgenes$SYMBOL,rownames(DFgenes))] 
      mygeneList <- egallgenes$stat
      names(mygeneList) <- egallgenes$ENTREZID
      mygeneList<-sort(mygeneList,decreasing = T)
      mygeneList <- mygeneList[!is.na(mygeneList)]
      head(mygeneList)
      length(mygeneList)
      
      # gseGO on all subontologies for differentially regulated genes
      for (m in c("BP", "MF", "CC")) {
        print(m)    
        GO.GSEA <- gseGO(geneList = mygeneList,
                         ont = m,
                         OrgDb = OrgDb,
                         keyType = "ENTREZID",
                         minGSSize = minGSSize,
                         maxGSSize = maxGSSize,
                         pvalueCutoff = pvalueCutoff,
                         pAdjustMethod = pAdjustMethod,
                         verbose = FALSE,
                         by= "fgsea")
        #ridgeplot first,grouped by gene set, density plots are generated by using the frequency of fold change values per gene within each set
        if (length(GO.GSEA@result$ID) >= 1) {
          rdgplot <- ridgeplot(GO.GSEA,showCategory = 20, label_format = 50)+ labs(x="enrichment distribution")
          # rdgplot <- rdgplot + xlim(NA,70)
          print(rdgplot)
          ggsave(filename = paste0(DEGFolder, "/GO_GSEA_ridgeplot_",m,".svg"), width = 16, height = 10)
          dev.off()
        }
        openxlsx::write.xlsx(GO.GSEA, file = paste0(DEGFolder,"/GO_GSEA",m,".xlsx"))
        GO.GSEA <- setReadable(GO.GSEA, OrgDb = OrgDb, keyType = "ENTREZID")
        openxlsx::write.xlsx(GO.GSEA, file = paste0(DEGFolder,"/GO_GSEA_readable",m,".xlsx"))
        GSEA.result[[species]][[paste0(marker.contrast,"_",m,"_result")]] <- GO.GSEA@result
        GOGSEA_File.Data.Frame <- as.data.frame(GO.GSEA)
        if (length(rownames(GOGSEA_File.Data.Frame)) < 1) {  ### if data frame is empty
          print("No significant terms detected")
        } else {
          (GOGSEA_File.Data.Frame[1:20,1:8])
          word.split <- strsplit(GO.GSEA@result$Description," ")
          
          ### linebreaks for plot, so the name doesnt become too long for the individual terms, a little bit outdated
          for (b in 1:length(word.split)) 
          {
            for (d in 1:length(word.split[[b]])) { ### checking for hyphen in words making them too long, erase them
              if (length(unlist(strsplit(word.split[[b]], " "))) < length(unlist(strsplit(word.split[[b]], "-"))))
              {
                nohyphen <- unlist(strsplit(word.split[[b]], "-"))
                word.split[[b]] <- nohyphen
              }
            }
            
            for (c in 1:length(word.split[[b]])) {  ### checking word length
              if (length(unlist(strsplit(word.split[[b]][c], split =""))) > 10)
              {
                word.split2 <- unlist(strsplit(word.split[[b]][c], split =""))
                word.split2 <- word.split2[1:10]
                word.split2 <- paste(word.split2, collapse = "")
                word.split[[b]][c] <- word.split2
              }
            }
            
            if (length(word.split[[b]]) > 4)
            {
              linebreaks <- word.split[[b]][1:4]
              linebreaks <- paste(linebreaks, collapse = " ")
              linebreaks2 <- word.split[[b]][5:length(word.split[[b]])]
              if (length(linebreaks2) > 4) 
              { 
                linebreaks2 <- linebreaks2[1:4]
                
                
              }
              linebreaks2 <- paste(linebreaks2, collapse = " ")
              word.split[[b]] <- paste(linebreaks,linebreaks2, sep = "\n" )
              
            } else (paste(word.split[[b]], collapse = " ") -> word.split[[b]])
            
          }
          GO.GSEA@result$Description <- word.split
          
          GO_GSEA <- GO.GSEA
          GOGSEA_File.Data.Frame <- as.data.frame(GO_GSEA)
          
          if (length(GO.GSEA@result$ID) > 20) {
            Set_Size <- GOGSEA_File.Data.Frame[1:20,]
          } else (Set_Size <- GOGSEA_File.Data.Frame[1:length(GO.GSEA@result$ID),])
          #bubbleplot
          pdf(paste0(DEGFolder, "/GO_GSEA_",m,".pdf"))
          ggplot(Set_Size, aes(x = NES , y = fct_reorder(unlist(Description), NES ))) +
            geom_point(aes(size = setSize, color = p.adjust)) +
            labs(size = "SetSize") +
            theme_bw(base_size = 14) +
            scale_colour_gradient(limits=c(min(GOGSEA_File.Data.Frame$p.adjust), max(GOGSEA_File.Data.Frame$p.adjust)), low="blue") +
            ylab(NULL) +
            expand_limits(x = c(min(Set_Size$NES)-(min(Set_Size$NES)/100*10), max(Set_Size$NES)+(max(Set_Size$NES)/100*10))) +
            ggtitle(paste0(m))-> plot.GSEA
          plot.GSEA + guides(size = guide_legend(order = 1)) -> plot.GSEA
          
          print(plot.GSEA)        
          dev.off()
          print(plot.GSEA)
          ggsave(filename = paste0(DEGFolder, "/GO_GSEA_",m,".svg"), width = 10, height = 10)
          # barplot
          GO_File <- Set_Size
          GO_File.Data.Frame <- as.data.frame(GO_File)
          
          if (length(rownames(GO_File.Data.Frame)) < 1) {  ### if data frame is empty
            warning("No significant terms detected")
          } else {
            GO_File.Data.Frame$NES <- as.numeric(GO_File.Data.Frame$NES)
            
            GO_File.Data.Frame2 <- GO_File.Data.Frame
            GO_File.Data.Frame2 <- GO_File.Data.Frame2[order(GO_File.Data.Frame2$p.adjust, decreasing = F),]
            GO_File.Data.Frame <-GO_File.Data.Frame2
            
            GO_File.Data.Frame$Description <- factor(GO_File.Data.Frame$Description, levels = unique(as.character(GO_File.Data.Frame$Description)))
            GO_File.Data.Frame <- transform(GO_File.Data.Frame, Description = reorder(Description, -p.adjust))
            if (length(rownames(GO_File.Data.Frame)) < 20) {k =length(rownames(GO_File.Data.Frame))} else {k=20}
            pdf(paste0(DEGFolder, "/GO_Barplot_",m,".pdf"))
            ggplot(GO_File.Data.Frame[1:k,], aes(reorder(GO_File.Data.Frame, p.adjust), x = Description , y = -log10(p.adjust))) +
              coord_flip() +
              geom_bar(aes( fill = setSize), stat="identity") +
              labs(size = "setSize") +
              theme_bw(base_size = 14) +
              theme(axis.text.y=element_text(size = 16)) +
              scale_colour_gradient(limits=c(min(GO_File.Data.Frame$setSize), max(GO_File.Data.Frame$setSize)), low="blue") +
              ylab("-log10(p.adjusted)") -> x
            x + guides(size = guide_legend(order = 1)) -> x
            print(x)
            dev.off()
            print(x)
            ggsave(filename = paste0(DEGFolder, "/GO_GSEA_Barplot_",m,".svg"), width = 10, height = 10)
            # gene concept plot
            if (length(rownames(GO_File.Data.Frame)) >= 5) {
              pdf(paste0(DEGFolder, "/GO_GSEA_GeneConcept_",m,".pdf"))          
              cnet.plot <- cnetplot(GO.GSEA, foldChange=mygeneList, node_label="all",cex_label_category = 1 , cex_label_gene = 0.6)
              print(cnet.plot)
              dev.off()
              print(cnet.plot)
              ggsave(filename = paste0(DEGFolder, "/GO_GSEA_GeneConcept_",m,".svg"), width = 16, height = 10)
              # circular gene concept plot
              pdf(paste0(DEGFolder, "/GO_GSEA_GeneConcept_circular_",m,".pdf")) 
              cnet.plot.circ <- cnetplot(GO.GSEA, foldChange = mygeneList, circular = TRUE, colorEdge = TRUE,cex_label_category = 1 , cex_label_gene = 0.6, showCategory = 4)
              print(cnet.plot.circ)
              dev.off()
              print(cnet.plot.circ)
              ggsave(filename = paste0(DEGFolder, "/GO_GSEA_GeneConcept_circular_",m,".svg"), width = 16, height = 10)
              
              
              
            }
          }
        }  
      }
    }  
  }
  return(GSEA.result) 
}

# Automated Annotation of Clusters for Seuratobjects 

#' @author Mariano Ruz Jurado
#' @param SeuratObject SeuratObject to annotate
#' @param ReferenceObject Annotated SeuratObject as Reference, yes it needs to be the same species...

DO.Annotation <- function(SeuratObject,ReferenceObject){
  require(SingleR)
  require(scran)
  require(scater)
  require(ggplotify)  
  #singleR works with SingleCellExperiment like objects, Seurat does have a function for it, also remove possible NAs from Reference
  ref <- as.SingleCellExperiment(ReferenceObject)
  ref <- ref[,!is.na(ref$cell_type)]
  
  SCE.SeuratObject <- as.SingleCellExperiment(SeuratObject)
  
  #Annotation with SingleR
  Cluster.Annotation <- SingleR(test=SCE.SeuratObject, ref=ref, labels=ref$cell_type, de.method="wilcox", clusters = SCE.SeuratObject$seurat_clusters)
  
  #Switch seurat cluster numbers with cell types
  SCE.SeuratObject$cell_type <- factor(
    SCE.SeuratObject$seurat_clusters,
    levels = rownames(Cluster.Annotation),
    labels = Cluster.Annotation$labels)
  
  SeuratObject$cell_type <- SCE.SeuratObject$cell_type
  Idents(SeuratObject) <- "cell_type"
  
  collected <- list()
  all.markers <- metadata(Cluster.Annotation)$de.genes
  empirical.markers <- findMarkers(SCE.SeuratObject, SCE.SeuratObject$cell_type, direction="up")
  
  
  for (i in unique(levels(SeuratObject$cell_type))) {
    
    markers <- unique(unlist(all.markers[[i]]))
    m <- match(markers, rownames(empirical.markers[[i]]))
    m <- markers[rank(m) <= 20] #top 20 upregulated marker genes found in both objects for celltype
    collected[[i]] <- as.grob(plotHeatmap(SCE.SeuratObject, order_columns_by="cell_type", features=m, main =i))
  }
  do.call(ggpubr::ggarrange, c(collected,ncol = 1))
  
  # Dim plots with annotated data
  p1<-DimPlot(SeuratObject, label = F, group.by = "species")
  p2<-DimPlot(SeuratObject, label = T)
  p3<-DimPlot(SeuratObject, label = F, group.by = "orig.ident")
  names(table(SeuratObject$cell_type))
  ggarrange.plot <- ggpubr::ggarrange(p1,p2,p3)
  print(ggarrange.plot)
  return(SeuratObject)
}

###protein sequence matching alignment
#' @author Mariano Ruz Jurado
#' @param mGene # Gene which will be matched against its possible orthologues in replacement
#' @param replacement # vector containing possible orthologues found by biomaRt
#' @param OrthologueList_allHuman # global orthologuelist which contains all Human genes and is filled with mouse orthologues
protein.matching <- function(mGene,replacement,OrthologueList_allHuman){
  
  require(mygene)
  require(UniprotR)
  require(RecordLinkage)
  require(Biostrings)
  
  outh <- try(queryMany(mGene, scopes = "symbol", fields= c("entrezgene", "uniprot"),species ="human"),silent = T)
  outm <- try(queryMany(replacement, scopes = "symbol", fields= c("entrezgene", "uniprot"),species ="mouse"),silent = T) # use try for some genes it gives errors
  
  if (!("try-error" %in% class(outm)) && !("try-error" %in% class(outh))) {
    
    #set a TrEMBL prot number or swiss prot number if TREMBL not available
    uniprot.mice <- list()
    for (i in 1:nrow(outm)) {
      if (!is.null(outm[i,]$uniprot.TrEMBL)) { # check if column is avalaible for TREMBL id if not set NA
        uniprot.mice[i] <- outm[i,]$uniprot.TrEMBL
      } else
      {
        uniprot.mice[i] <- NA
      }
      if (is.na(uniprot.mice[i]) && !is.null(outm[i,]$uniprot.Swiss.Prot)) { # check if swiss prot number available
        uniprot.mice[i] <- outm[i,]$uniprot.Swiss.Prot
      }
    }
    
    #set not found entries with NA
    for (i in 1:length(uniprot.mice)) {
      if (!is.null(uniprot.mice[[i]][1]) && !is.na(uniprot.mice[[i]][1])) {
        uniprot.mice[[i]] <- uniprot.mice[[i]][1]
      } else{
        uniprot.mice[[i]] <- NA # set NA to not found in database
        replacement[i] <- NA # set NA to these Gene Symbols as well
        
      }
    }
    #delete NAs
    replacement <- replacement[!is.na(replacement)]
    uniprot.mice <- uniprot.mice[!is.na(uniprot.mice)]
    
    #get sequences based on uniprot IDs
    #human
    sequences.h <- suppressWarnings(GetSequences(unlist(outh$uniprot.Swiss.Prot), directorypath = NULL))
    human.sequences <- sequences.h$Sequence
    if (plyr::empty(sequences.h)) {
      sequences.h <- suppressWarnings(GetSequences(unlist(outh$uniprot.TrEMBL), directorypath = NULL))
      human.sequences <- sequences.h$Sequence
    }
    #mice
    sequences.m <- suppressWarnings(GetSequences(unlist(uniprot.mice), directorypath = NULL))
    mice.sequences <- sequences.m$Sequence
    
    #pairwise alignment for checking sequence similarity, returning a score which is then used for determining best orthologue
    if (!is.null(human.sequences) && !is.null(mice.sequences)) {
      local.Align.list <- list()
      for (k in 1:length(mice.sequences)) {
        localAlign <- pairwiseAlignment(human.sequences,mice.sequences[k], substitutionMatrix="BLOSUM50" , gapOpening = 0, gapExtension = 8)
        local.Align.list[[k]] <- localAlign@score
      }
      names(local.Align.list) <- replacement #set replacement names
      local.Align.list <- local.Align.list[order(-unlist(local.Align.list))] #order by score
      replacement.hit <- names(local.Align.list[1]) # first entry now has the highest score
      j<-1 # set counter for while
      while (replacement.hit %in% OrthologueList_allHuman$MouseGene && !is.na(replacement.hit)) { # check if already in DF, if yes then take second hit
        replacement.hit <- names(local.Align.list[j+1])
        j<-j+1
      }
      OrthologueList_allHuman[OrthologueList_allHuman$HGNC.symbol==mGene,]$MouseGene=replacement.hit # set orthologue
      
    }
  }
  else
  {
    human.sequences <- NA
    mice.sequences <- NA
  }
  return(list(OrthologueList_allHuman,human.sequences,mice.sequences))
}



###nucleotide sequence comparison alignment
#' @author Mariano Ruz Jurado
#' @param mGene # Gene which will be matched against its possible orthologues in replacement
#' @param replacement # vector containing possible orthologues found by biomaRt
#' @param OrthologueList_allHuman # global orthologuelist which contains all Human genes and is filled with mouse orthologues
nucleotide.matching <- function(mGene,replacement,OrthologueList_allHuman){
  
  require(rentrez)
  require(mygene)
  require(stringr)
  require(Biostrings)
  out <- try(queryMany(mGene, scopes = "symbol", fields= c("entrezgene", "uniprot"),species ="human"),silent = T)
  outm <- try(queryMany(replacement, scopes = "symbol", fields= c("entrezgene", "uniprot"),species ="mouse"),silent = T) # use try for some genes it gives errors
  
  if (!("try-error" %in% class(outm)) && !("try-error" %in% class(out))) {
    #human
    linked_seq_ids <- entrez_link(dbfrom = "gene", id=out$entrezgene, db="nuccore")
    linked_transcripts <- linked_seq_ids$links$gene_nuccore_refseqrna
    head(linked_transcripts)
    if (!is.null(linked_transcripts)) { #check if it hits a sequence human
      
      all_recs <- entrez_fetch(db="nuccore", id=linked_transcripts, rettype = "fasta")
      sequences <- unlist(strsplit(all_recs, split=">"))
      sequences.orth <- sequences[grep("NM",sequences)] # NM is for non predicted mRNAs
      if (identical(sequences.orth, character(0)) == T) {
        sequences.orth <- sequences[grep("NR", sequences)] ## check if it is non coding RNA if mRNA is 0
      }
      if (length(sequences.orth) > 1) {
        sequences.orth <- sequences.orth[grep("transcript variant 1", sequences.orth)]
      }
      
      
      sequences.orth <- unlist(strsplit(sequences.orth, split="\n"))
      Non.variant.h <- sequences.orth[2:(length(sequences.orth)-1)]
      Non.variant.h <- paste0(Non.variant.h, collapse="")
    }
    
    #mouse
    if (!is.null(outm$entrezgene)) {
      linked_seq_ids <- entrez_link(dbfrom = "gene", id=outm$entrezgene, db="nuccore")
      linked_transcripts <- linked_seq_ids$links$gene_nuccore_refseqrna
      head(linked_transcripts)
    }
    if (!is.null(linked_transcripts)) { #check if it hits a sequence mice
      
      all_recs <- entrez_fetch(db="nuccore", id=linked_transcripts, rettype = "fasta")
      sequences <- unlist(strsplit(all_recs, split=">"))
      sequences.orth <- sequences[grep("NM",sequences)] # NM is for non predicted mRNAs
      if (identical(sequences.orth, character(0)) == T) {
        sequences.orth <- sequences[grep("NR", sequences)] ## check if there is non coding RNA if mRNA is 0
      }
      if (identical(sequences.orth, character(0)) == T) {
        sequences.orth <- sequences[grep("XM", sequences)] ## check if there is a predicted coding mRNA
      }
      if (identical(sequences.orth, character(0)) == T) {
        sequences.orth <- sequences[grep("XR", sequences)] ## check if there is a predicted non coding mRNA
      }
      
      orthologues.sequences <- as.list(sequences.orth)
      for (i in 1:length(orthologues.sequences)) {
        orthologues.seq <- unlist(strsplit(orthologues.sequences[[i]], split="\n"))
        orthologues.seq.2 <- orthologues.seq[2:length(orthologues.seq)]
        orthologues.sequences[[i]] <- paste0(orthologues.seq.2, collapse="")
        names(orthologues.sequences)[i] <- orthologues.seq[1]
      }
    }
    if (!is.null(linked_transcripts)) { # no transcript sequences, no alignment
      #alignment
      local.Align.list <- list()
      mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = FALSE, type = "DNA")
      for (k in 1:length(orthologues.sequences)) {
        localAlign <- pairwiseAlignment(
          AAString(Non.variant.h),
          AAString(orthologues.sequences[[k]]),
          type="local",
          substitutionMatrix=mat , gapOpening = 5, gapExtension = 2)
        local.Align.list[[k]] <- localAlign@score
        #get gene name from string, stands between "()"
        names(local.Align.list)[k] <- str_match(names(orthologues.sequences)[k],"[(](.*?)[)]")[,2]
      }
      local.Align.list <- local.Align.list[order(-unlist(local.Align.list))] #order by score
      replacement.hit <- names(local.Align.list[1]) # first entry now has the highest score
      j<-1 # set counter for while
      while (replacement.hit %in% OrthologueList_allHuman$MouseGene && !is.na(replacement.hit)) { # check if already in DF, if yes then take second hit
        replacement.hit <- names(local.Align.list[j+1])
        j<-j+1
      }
      OrthologueList_allHuman[OrthologueList_allHuman$HGNC.symbol==mGene,]$MouseGene=replacement.hit # set orthologue
    }
  }
  else
  {
    linked_transcripts <- NA
  }
  return(list(OrthologueList_allHuman, linked_transcripts))
}
