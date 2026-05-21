# tutorial code adapted from link below:
# https://www.archrproject.com/articles/Articles/tutorial.html

library(ArchR)
library(parallel)
#library(BiocManager)

## Setting ArchRLocking to TRUE.
addArchRLocking(locking = TRUE)

set.seed(1)

#Gets tutorial data.
inputFiles <- getTutorialData("Hematopoiesis")

# If no internet, run the following instead to get names of input files
#inputFiles <- c("HemeFragments/scATAC_BMMC_R1.fragments.tsv.gz",
 #               "HemeFragments/scATAC_CD34_BMMC_R1.fragments.tsv.gz",
  #              "HemeFragments/scATAC_PBMC_R1.fragments.tsv.gz" )

#adding names to files
#names(inputFiles) <- c("scATAC_BMMC_R1", "scATAC_CD34_BMMC_R1", "scATAC_PBMC_R1")

head(inputFiles)

#Create directory folder to download database R package
#dir.create("/home/lionavu/projects/def-itobias/BINF_6999/R_libs", recursive = TRUE)

#Install package directly from folder.
# first ran in terminal the following
# wget https://bioconductor.org/packages/release/data/annotation/src/contrib/BSgenome.Hsapiens.UCSC.hg19_1.4.3.tar.gz

#install.packages("/home/lionavu/projects/def-itobias/BINF_6999/BSgenome.Hsapiens.UCSC.hg19_1.4.3.tar.gz",
#               lib = "/home/lionavu/projects/def-itobias/BINF_6999/R_libs",
#               repos = NULL,
#               type = "source")

#Load library from folder
library(BSgenome.Hsapiens.UCSC.hg19,
        lib.loc = "/home/lionavu/projects/def-itobias/BINF_6999/R_libs")

## Setting default genome to Hg19, for some reason, it did not work without downloading the file above
addArchRGenome("hg19")

#Set parallel threads to capability to 16
addArchRThreads(threads = 16)

#Create arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE)

# Line of code here prevents parallelization issues by setting thread to 1
addArchRThreads(threads = 1)

#Identify doublets in data, following had to be run on the rstudio server to work, HPC did not work...
doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1)

#Create arrow project
projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "HemeTutorial",
  copyArrows = TRUE)

#save the following to transfer file to HPC
saveRDS(projHeme1, file = "projHeme1.RDS", compress = FALSE)

#Run the following line on HPC
projHeme1 <- readRDS("projHeme1.RDS")

#Look at contents of ArchR project
projHeme1

#look at memory size of project
paste0("Memory Size = ", round(object.size(projHeme1) / 10^6, 3), " MB")
#[1] "Memory Size = 37.477 MB"

#Check which data matrices are available in ArchR project
getAvailableMatrices(projHeme1)

#Look at all metadata in cellColData available 
head(projHeme1@cellColData)

#Access cellColData by using the $ or the @ operator
head(projHeme1@cellNames)
head(projHeme1@cellColData@rownames)

#Access sample names by using the $ operator
head(projHeme1$Sample)

# Access the TSS Enrichment Scores for each cell
head(projHeme1$TSSEnrichment)
quantile(projHeme1$TSSEnrichment)

getArrowFiles(projHeme1)
#directory of arrow files are outdated and is from the R server, must change to match to HPC directories

#Change arrow file directory in metadata to match with directory on HPC cluster
projHeme1@sampleColData$ArrowFiles[1] <- "/home/lionavu/projects/def-itobias/BINF_6999/archr_tutorial/scATAC_BMMC_R1.arrow" 
projHeme1@sampleColData$ArrowFiles[2] <- "/home/lionavu/projects/def-itobias/BINF_6999/archr_tutorial/scATAC_CD34_BMMC_R1.arrow"
projHeme1@sampleColData$ArrowFiles[3] <- "/home/lionavu/projects/def-itobias/BINF_6999/archr_tutorial/scATAC_PBMC_R1.arrow"

#create directory to copy archr subsetted project into, else you will get a recursive error
#dir.create("/home/lionavu/projects/def-itobias/BINF_6999/archr_tutorial/ArchR_subset_1", recursive = TRUE)

#Change output file directory in metadata to match with directory on HPC cluster
projHeme1@projectMetadata@listData$outputDirectory <- "/home/lionavu/projects/def-itobias/BINF_6999/archr_tutorial/ArchR_subset_1"

#create archr subset
projSubset <- subsetArchRProject(
  ArchRProj = projHeme1,
  cells = projHeme1$cellNames[idxSample],
  outputDirectory = "ArchRSubset",
  dropCells = TRUE,
  force = TRUE)

#subset by first 100 cells
projHeme1[1:100, ]

projHeme1[projHeme1$cellNames[1:100], ]

#plot ridge plots
p2 <- plotGroups(
    ArchRProj = projHeme1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
    baseSize = 10,
  addBoxPlot = TRUE,)

# ggsave(p1, filename = "figs/ridge.png")

p2 <- plotGroups(
    ArchRProj = projHeme1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
    baseSize = 10,
  addBoxPlot = TRUE,)
# ggsave(p2, filename = "figs/violin_plots.png")

#Example 3. Make a ridge plot for each sample for the log10(unique nuclear fragments).

p3 <- plotGroups(
    ArchRProj = projHeme1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges",
    baseSize = 10)
## 1

p3
## Picking joint bandwidth of 0.0506

#Example 4. Make a violin plot for each sample for the log10(unique nuclear fragments).
p4 <- plotGroups(
    ArchRProj = projHeme1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
    alpha = 0.4,
    baseSize = 10,
  addBoxPlot = TRUE)
## 1

p4

#Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles.
p5 <- plotFragmentSizes(ArchRProj = projHeme1)
p6 <- plotTSSEnrichment(ArchRProj = projHeme1)

#Saving arch R project
projHeme1 <- saveArchRProject(ArchRProj = projHeme1, outputDirectory = "Save-ProjHeme1", load = TRUE)

#loading the project
projHeme1 <- loadArchRProject(path = "./Save-ProjHeme1")

#filtering for doublets
projHeme2 <- filterDoublets(projHeme1)
#Filtering 410 cells from ArchRProject!
	#scATAC_BMMC_R1 : 243 of 4932 (4.9%)
	#scATAC_CD34_BMMC_R1 : 107 of 3275 (3.3%)
	#scATAC_PBMC_R1 : 60 of 2453 (2.4%)

#compare number of cells pre and post doublet filtering
length(getCellNames(ArchRProj = projHeme1))
#[1] 10660
length(getCellNames(ArchRProj = projHeme2))
#[1] 10250

#perform iterative LSI for dimensionality reduction, ensure threads = 1!!!!! HPC with apptainer again did not work, but worked on Rserver
#may be an issue with fortran which irbla relies on and is breaking the lsi computation along with the open blast system that apptainer runs on
#setting this to 1 thread in R may help
#RhpcBLASctl::blas_set_num_threads(1)
#RhpcBLASctl::omp_set_num_threads(1)

projHeme2 <- addIterativeLSI(
    ArchRProj = projHeme2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)

#Batch effect correction with Harmony
projHeme2 <- addHarmony(
    ArchRProj = projHeme2,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample")

#Add clustering using Seurat
projHeme2 <- addClusters(input = projHeme2,
                         reducedDims = "IterativeLSI",
                         method = "Seurat",
                         name = "Clusters",
                         resolution = 0.8)

#Look at clusters using the $ accessor which shows the cluster ID for each single cell.
head(projHeme2$Clusters)
# [1] "C4" "C6" "C9" "C9" "C9" "C4"

# Tabulate the number of cells present in each cluster:
table(projHeme2$Clusters)
#  C1  C10  C11  C12   C2   C3   C4   C5   C6   C7   C8   C9 
#1495  422  312  406 1118  909 1089 1362 1240  886  656  355 

#Create confusion matrix to see which cells got put in which cluster
cM <- confusionMatrix(paste0(projHeme2$Clusters), paste0(projHeme2$Sample))
cM
#12 x 3 sparse Matrix of class "dgCMatrix"
#    scATAC_BMMC_R1 scATAC_CD34_BMMC_R1 scATAC_PBMC_R1
#C4             293                 795              1
#C6            1201                   .             39
#C9             355                   .              .
#C10            254                   5            163
#C1            1459                   4             32
#C7             305                   .            581
#C3             221                 678             10
#C8              98                   1            557
#C11            161                 141             10
#C2             117                   1           1000
#C5             133                1229              .
#C12             92                 314              .

#create heatmap from confusion matrix
library(pheatmap)

cM <- cM / Matrix::rowSums(cM)
p_heatmap <- pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black")

# ggsave(p_heatmap, file = "Save-ProjHeme3/Plots/heatmap_confusion_mat.png")

#ArchR supports clustering using scran, install these packages first:
#edgeR v4.2.2, locfit, metapod, assorthead, 

#projHeme2 <- addClusters(
 #   input = projHeme2,
 #   reducedDims = "IterativeLSI",
 #   method = "scran",
 #   name = "ScranClusters",
 #   k = 15)

#Generate UMAP
projHeme2 <- addUMAP(
    ArchRProj = projHeme2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine")

#plot UMAP using plotEmbedding function
p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

#save plot above 
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

projHeme2 <- addTSNE(
    ArchRProj = projHeme2, 
    reducedDims = "IterativeLSI", 
    name = "TSNE", 
    perplexity = 30)

p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")

#save tSNE plot
plotPDF(p1,p2, name = "Plot-TSNE-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

#repeat the above but with batch correction for UMAP
projHeme2 <- addUMAP(
    ArchRProj = projHeme2, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

#identify marker features
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

#visualize markers
markerGenes  <- c(
    "CD34", #Early Progenitor
    "GATA1", #Erythroid
    "PAX5", "MS4A1", "EBF1", "MME", #B-Cell Trajectory
    "CD14", "CEBPB", "MPO", #Monocytes
    "IRF8", 
    "CD3D", "CD8A", "TBX21", "IL7R" #TCells )
	
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE)

#plot heatmap
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
heatmapGS@row_order <- c(10,11,12,3,4,5,1,2,6,8,7,9)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme2, addDOC = FALSE)

markerGenes  <- c(
    "CD34",  #Early Progenitor
    "GATA1", #Erythroid
    "PAX5", "MS4A1", "MME", #B-Cell Trajectory
    "CD14", "MPO", #Monocytes
    "CD3D", "CD8A"#TCells
  )

p <- plotEmbedding(
    ArchRProj = projHeme2, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL)

p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

plotPDF(plotList = p, 
    name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf", 
    ArchRProj = projHeme2, 
    addDOC = FALSE, width = 5, height = 5)

p <- plotEmbedding(
    ArchRProj = projHeme2, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHeme2))

#Rearrange for grid plotting
p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

#Add module scores
features <- list(
  BScore = c("MS4A1", "CD79A", "CD74"),
  TScore = c("CD3D", "CD8A", "GZMB", "CCR7", "LEF1"))
	
projHeme2 <- addModuleScore(projHeme2,
    useMatrix = "GeneScoreMatrix",
    name = "Module",
    features = features)

p5 <- plotEmbedding(projHeme2,
    embedding = "UMAP",
    colorBy = "cellColData",
    name="Module.BScore",
    imputeWeights = getImputeWeights(projHeme2))

p6 <- plotEmbedding(projHeme2,
    embedding = "UMAP",
    colorBy = "cellColData",
    name="Module.TScore",
    imputeWeights = getImputeWeights(projHeme2))

plotPDF(ggAlignPlots(p5,p6,draw=F,type="h"))

#Track Plotting with ArchRBrowser
set.seed(1)
plot_browser <- plotBrowserTrack(
    ArchRProj = projHeme2, 
    groupBy = "Clusters", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000)

#plot the actual track 
grid::grid.newpage()
grid::grid.draw(p$CD14)

#save a multi-page PDF with a single page for each gene locus 
plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes.pdf", 
    ArchRProj = projHeme2, 
    addDOC = FALSE, width = 5, height = 5)

#Defining clusters with scRNA seq data
#run wget https://jeffgranja.s3.amazonaws.com/ArchR/TestData/scRNA-Hematopoiesis-Granja-2019.rds on terminal beforehand
	
seRNA <- readRDS("scRNA-Hematopoiesis-Granja-2019.rds")
seRNA

#eda of scrna seq data
colnames(colData(seRNA))
table(colData(seRNA)$BioClassification)
projHeme2 <- loadArchRProject(path = "./Save-ProjHeme1-step-10")
seRNA <- readRDS("scRNA-Hematopoiesis-Granja-2019.rds")
seRNA

#generating an unconstrained integration, works on Rstudio server with docker image but not in HPC
projHeme2 <- addGeneIntegrationMatrix(
    ArchRProj = projHeme2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = FALSE,
    groupRNA = "BioClassification",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un")

#check if gene score matrix was added
getAvailableMatrices(projHeme2)

#performing constrained integration
cM <- as.matrix(confusionMatrix(projHeme2$Clusters, projHeme2$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments

# look at the cell type labels from our scRNA-seq data 
unique(unique(projHeme2$predictedGroup_Un))

#From scRNA
cTNK <- paste0(paste0(19:25), collapse="|")
cTNK

cNonTNK <- paste0(c(paste0("0", 1:9), 10:13, 15:18), collapse="|")
cNonTNK

# For T cells and NK cells clusters
clustTNK <- rownames(cM)[grep(cTNK, preClust)]
clustTNK

#identifies non T and non NK cells
clustNonTNK <- rownames(cM)[grep(cNonTNK, preClust)]
clustNonTNK

#identify T and NK cells in RNA seq data
rnaTNK <- colnames(seRNA)[grep(cTNK, colData(seRNA)$BioClassification)]
head(rnaTNK)

#identify the Non-T cell Non-NK cell
rnaNonTNK <- colnames(seRNA)[grep(cNonTNK, colData(seRNA)$BioClassification)]
head(rnaNonTNK)

#nested list for constrained integration
groupList <- SimpleList(
    TNK = SimpleList(
        ATAC = projHeme2$cellNames[projHeme2$Clusters %in% clustTNK],
        RNA = rnaTNK),
    NonTNK = SimpleList(
        ATAC = projHeme2$cellNames[projHeme2$Clusters %in% clustNonTNK],
        RNA = rnaNonTNK))

#perform constrained integration
projHeme2 <- addGeneIntegrationMatrix(
    ArchRProj = projHeme2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = FALSE, 
    groupList = groupList,
    groupRNA = "BioClassification",
    nameCell = "predictedCell_Co",
    nameGroup = "predictedGroup_Co",
    nameScore = "predictedScore_Co")

#compare both constrained and unconstrained integration
pal <- paletteDiscrete(values = colData(seRNA)$BioClassification)
pal

#plot 
p1 <- plotEmbedding(
    projHeme2, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un", 
    pal = pal)

p2 <- plotEmbedding(
    projHeme2, 
    colorBy = "cellColData", 
    name = "predictedGroup_Co", 
    pal = pal)

#save plot
plotPDF(p1,p2, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

#projHeme2 <- saveArchRProject(ArchRProj = projHeme2, outputDirectory = "Save-ProjHeme2-step-10-1", load = TRUE)
#projHeme2 <- loadArchRProject(path = "./Save-ProjHeme2-step-10-1")

#Adding Pseudo-scRNA-seq profiles for each scATAC-seq cell, this can be computationally demanding! works on server but not on HPC apptainer
projHeme3 <- addGeneIntegrationMatrix(
    ArchRProj = projHeme2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = TRUE,
    force= TRUE,
    groupList = groupList,
    groupRNA = "BioClassification",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore")

getAvailableMatrices(projHeme3)

saveArchRProject(ArchRProj = projHeme3, outputDirectory = "Save-ProjHeme2-step-10-1.1", load = TRUE)

#added impute weights to our project
projHeme3 <- addImputeWeights(projHeme3)

#gene expression UMAP overlay, defining cell types
markerGenes  <- c(
    "CD34", #Early Progenitor
    "GATA1", #Erythroid
    "PAX5", "MS4A1", #B-Cell Trajectory
    "CD14", #Monocytes
    "CD3D", "CD8A", "TBX21", "IL7R") #TCells

p1 <- plotEmbedding(
    ArchRProj = projHeme3, 
    colorBy = "GeneIntegrationMatrix", 
    name = markerGenes, 
    continuousSet = "horizonExtra",
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHeme3))

#now with the gene score
p2 <- plotEmbedding(
    ArchRProj = projHeme3, 
    colorBy = "GeneScoreMatrix", 
    continuousSet = "horizonExtra",
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHeme3))

#plot all marker genes
p1c <- lapply(p1, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})

p2c <- lapply(p2, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
	
#show all plots for gene score and gene expression
do.call(cowplot::plot_grid, c(list(ncol = 3), p1c))
do.call(cowplot::plot_grid, c(list(ncol = 3), p2c))

#save plots
plotPDF(plotList = p1, 
    name = "Plot-UMAP-Marker-Genes-RNA-W-Imputation.pdf", 
    ArchRProj = projHeme3, 
    addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = projHeme3, outputDirectory = "./Save-ProjHeme2-step-10-2", load = TRUE)

cM <- confusionMatrix(projHeme3$Clusters, projHeme3$predictedGroup)
labelOld <- rownames(cM)

labelNew <- colnames(cM)[apply(cM, 1, which.max)]

remapClust <- c(
    "01_HSC" = "Progenitor",
    "02_Early.Eryth" = "Erythroid",
    "03_Late.Eryth" = "Erythroid",
    "04_Early.Baso" = "Basophil",
    "05_CMP.LMPP" = "Progenitor",
    "06_CLP.1" = "CLP",
    "07_GMP" = "GMP",
    "08_GMP.Neut" = "GMP",
    "09_pDC" = "pDC",
    "10_cDC" = "cDC",
    "11_CD14.Mono.1" = "Mono",
    "12_CD14.Mono.2" = "Mono",
    "13_CD16.Mono" = "Mono",
    "15_CLP.2" = "CLP",
    "16_Pre.B" = "PreB",
    "17_B" = "B",
    "18_Plasma" = "Plasma",
    "19_CD8.N" = "CD8.N",
    "20_CD4.N1" = "CD4.N",
    "21_CD4.N2" = "CD4.N",
    "22_CD4.M" = "CD4.M",
    "23_CD8.EM" = "CD8.EM",
    "24_CD8.CM" = "CD8.CM",
    "25_NK" = "NK")
	
remapClust <- remapClust[names(remapClust) %in% labelNew]

labelNew2 <- mapLabels(labelNew, oldLabels = names(remapClust), newLabels = remapClust)
projHeme3$Clusters2 <- mapLabels(projHeme3$Clusters, newLabels = labelNew2, oldLabels = labelOld)

p1 <- plotEmbedding(projHeme3, colorBy = "cellColData", name = "Clusters2")
projHeme3 <- saveArchRProject(ArchRProj = projHeme3, outputDirectory = "Save-ProjHeme3", load = TRUE)

#making pseudo bulk replicaatesm this takes a loong time
projHeme4 <- addGroupCoverages(ArchRProj = projHeme3, groupBy = "Clusters2")

groups <- addGroupCoverages(ArchRProj = projHeme3, groupBy = "Clusters2", returnGroups = TRUE)

#saveArchRProject(ArchRProj = projHeme4, outputDirectory = "Save-ProjHeme2-step-11-2", load = TRUE)
	
#projHeme4 <- loadArchRProject(path = "/.Save-ProjHeme2-step-11-2", load =TRUE)
							  
#now to calling peaks with macs2
#run pip install macs2 first on terminal!
#make sure to load the genome beforehand! else r will complain ...

#find macs directory
pathToMacs2 <- findMacs2()
#Searching For MACS2..
#Found with $PATH at /usr/local/bin/macs2
	
projHeme4 <- addReproduciblePeakSet(
    ArchRProj = projHeme4, 
    groupBy = "Clusters2", 
    pathToMacs2 = pathToMacs2)









	


