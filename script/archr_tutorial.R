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








