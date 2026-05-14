library(ArchR)
library(parallel)
#library(BiocManager)

## Setting ArchRLocking to TRUE.
addArchRLocking(locking = TRUE)

set.seed(1)

#inputFiles <- getTutorialData("Hematopoiesis")
inputFiles <- c("HemeFragments/scATAC_BMMC_R1.fragments.tsv.gz",
                "HemeFragments/scATAC_CD34_BMMC_R1.fragments.tsv.gz",
                "HemeFragments/scATAC_PBMC_R1.fragments.tsv.gz" )

names(inputFiles) <- c("scATAC_BMMC_R1", "scATAC_CD34_BMMC_R1", "scATAC_PBMC_R1")

head(inputFiles)
inputFiles_example <- c("/path/to/fragFile1.tsv.gz", "/path/to/fragFile2.tsv.gz")

#adding names to files
names(inputFiles_example) <- c("Sample1","Sample2")

#Create directory folder to download database R package
#dir.create("/home/lionavu/projects/def-itobias/BINF_6999/R_libs", recursive = TRUE)

#install.packages("/home/lionavu/projects/def-itobias/BINF_6999/BSgenome.Hsapiens.UCSC.hg19_1.4.3.tar.gz",
#               lib = "/home/lionavu/projects/def-itobias/BINF_6999/R_libs",
#               repos = NULL,
#               type = "source")

#Load library from folder
library(BSgenome.Hsapiens.UCSC.hg19,
        lib.loc = "/home/lionavu/projects/def-itobias/BINF_6999/R_libs")

## Setting default genome to Hg19
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

addArchRThreads(threads = 1)

#Identify doublets in data 
doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1)

