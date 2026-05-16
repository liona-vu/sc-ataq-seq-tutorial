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

#Install package directly from folder
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

#Change output file directory in metadata to match with directory on HPC cluster
projHeme1@projectMetadata@listData$outputDirectory <- "/home/lionavu/projects/def-itobias/BINF_6999/archr_tutorial/ArchR_subset_1"

#create directory to copy archr subsetted project into, else you will get a recursive error
#dir.create("/home/lionavu/projects/def-itobias/BINF_6999/archr_tutorial/ArchR_subset_1", recursive = TRUE)

#create archr subset
projSubset <- subsetArchRProject(
  ArchRProj = projHeme1,
  cells = projHeme1$cellNames[idxSample],
  outputDirectory = "ArchRSubset",
  dropCells = TRUE,
  force = TRUE)

projHeme1[1:100, ]

projHeme1[projHeme1$cellNames[1:100], ]


