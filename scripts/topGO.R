library(topGO)
library(stats)
library(knitr)


setwd('/Users/mmayers/projects/n15_mice/data/')

data.dir <- '../data'

gene2GOID <- readMappings('clusterID2GO.map')
unenr_pvals <- read.csv('unenriched_pvals.csv')
unenr_pvals$log2FoldChange <- log2(unenr_pvals$ratio)
unenr_pvals$padj <- p.adjust(unenr_pvals$p_value, method="BH")
GO2geneMF <- annFUN.gene2GO('MF', gene2GO=gene2GOID)
GO2geneBP <- annFUN.gene2GO('BP', gene2GO=gene2GOID)

# do some stuff with fold changes
unenr_pvals$adjFC <- unenr_pvals$log2FoldChange
n <- sum(unenr_pvals$padj > .05 & abs(unenr_pvals$log2FoldChange) > 2)
set.seed(2222)
unenr_pvals[unenr_pvals$padj > .05 & abs(unenr_pvals$log2FoldChange) > 2, ]$adjFC <- runif(n, -1, 1)
hist(unenr_pvals$log2FoldChange, 15)
hist(unenr_pvals$adjFC, 15)

# Metadata on Cluster Info - Use this for filtering of certain types of proteins
sampleTable <- read.csv(file.path(data.dir, 'filt_metadata.csv'), row.names=1)
sampleTable$technical <- as.logical(sampleTable$technical)
locusTable <- read.csv(file.path(data.dir, 'loci_annot.csv'))


# Look at Bacterial Proteins only (Mouse Human subset seems to be a bust, can't figure out
# good way to account for tiny background)
tmp <- locusTable$mouse_human == 'False' & locusTable$lca != '35823' # Filtering out anthrospira differences 35823
tmp[is.na(tmp)] <- TRUE
bac.prots <- locusTable[tmp, ]$X
bac.pvals <- unenr_pvals[unenr_pvals$X %in% bac.prots, ]
geneList.bac <- bac.pvals$padj
geneList.bac <- bac.pvals$adjFC
names(geneList.bac) <- bac.pvals$X

topDiffGenes <- function(score) {return(score >= 2)} # scoring function, anything with > 2 logfc is sig anyway

## Look at what's inside Tcell vs RAG
geneList.bac <- -1*geneList.bac

ontology <- 'MF'
GOdata.bac <- new('topGOdata', ontology = ontology, allGenes = geneList.bac, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.bac.ks <- runTest(GOdata.bac, statistic = 'ks')
GenTable(GOdata.bac, result.bac.ks, topNodes = 10)
showGroupDensity(GOdata.bac, "GO:0008184", ranks=TRUE) #Glycogen phosphorlase
showGroupDensity(GOdata.bac, "GO:0004618", ranks=TRUE) #phosphoglycerate kinase

ontology <- 'BP'
GOdata.bac <- new('topGOdata', ontology = ontology, allGenes = geneList.bac, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.bac.ks <- runTest(GOdata.bac, statistic = 'ks')
GenTable(GOdata.bac, result.bac.ks, topNodes = 10)
showGroupDensity(GOdata.bac, "GO:0006096", ranks=TRUE) # Glycolitc Process
showGroupDensity(GOdata.bac, "GO:0005975", ranks=TRUE) # Carbohydydrate metab
# Mostly Carb metabolism?  lipid biosynthesis?Nah, those look crappy

# RAG vs Tcell
geneList.bac <- -1*geneList.bac

ontology <- 'MF'
GOdata.bac <- new('topGOdata', ontology = ontology, allGenes = geneList.bac, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.bac.ks <- runTest(GOdata.bac, statistic = 'ks')
GenTable(GOdata.bac, result.bac.ks, topNodes = 10)

ontology <- 'BP'
GOdata.bac <- new('topGOdata', ontology = ontology, allGenes = geneList.bac, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.bac.ks <- runTest(GOdata.bac, statistic = 'ks')
GenTable(GOdata.bac, result.bac.ks, topNodes = 10)
showGroupDensity(GOdata.bac, "GO:0005978", ranks=TRUE) #Glycogen biosynthetic process


#BioGlyCMK Rag vs Tcell
enr_pvals <- read.csv(file.path(data.dir, 'enriched_pvals.csv'))
enr_pvals$log2FoldChange <- log2(enr_pvals$ratio)
enr_pvals$padj <- p.adjust(enr_pvals$p_value, method="BH")

# do some stuff with fold changes
enr_pvals$adjFC <- enr_pvals$log2FoldChange
n <- sum(enr_pvals$padj > .05 & abs(enr_pvals$log2FoldChange) > 2)
set.seed(1234)
enr_pvals[enr_pvals$padj > .05 & abs(enr_pvals$log2FoldChange) > 2, ]$adjFC <- runif(n, -1, 1)
hist(enr_pvals$log2FoldChange, 15)
hist(enr_pvals$adjFC, 15)

# Bacterial proteins only
bac.enr.pvals <- enr_pvals[enr_pvals$X %in% bac.prots, ]
#geneList.bac <- bac.pvals$padj
geneList.bac.enr <- bac.enr.pvals$adjFC
names(geneList.bac.enr) <- bac.enr.pvals$X

## Look at what's inside 

# Tcell vs RAG
geneList.bac.enr <- -1*geneList.bac.enr # KS goes from - to +
                                        # do this to put upreg in Tcell first
ontology <- 'MF'
GOdata.bac.enr <- new('topGOdata', ontology = ontology, allGenes = geneList.bac.enr, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.bac.enr.ks <- runTest(GOdata.bac.enr, statistic = 'ks')
GenTable(GOdata.bac.enr, result.bac.enr.ks, topNodes = 10)

ontology <- 'BP'
GOdata.bac.enr <- new('topGOdata', ontology = ontology, allGenes = geneList.bac.enr, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.bac.enr.ks <- runTest(GOdata.bac.enr, statistic = 'ks')
GenTable(GOdata.bac.enr, result.bac.enr.ks, topNodes = 10)
showGroupDensity(GOdata.bac.enr, "GO:0006508", ranks = TRUE) #proteolysis

# RAG vs TCell
ontology <- 'MF'
geneList.bac.enr <- -1*geneList.bac.enr 

GOdata.bac.enr <- new('topGOdata', ontology = ontology, allGenes = geneList.bac.enr, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.bac.enr.ks <- runTest(GOdata.bac.enr, statistic = 'ks')
GenTable(GOdata.bac.enr, result.bac.enr.ks, topNodes = 10)

ontology <- 'BP'
GOdata.bac.enr <- new('topGOdata', ontology = ontology, allGenes = geneList.bac.enr, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.bac.enr.ks <- runTest(GOdata.bac.enr, statistic = 'ks')
GenTable(GOdata.bac.enr, result.bac.enr.ks, topNodes = 10)
showGroupDensity(GOdata.bac.enr, "GO:0006091", ranks=TRUE) #generation of precursor metabolites
showGroupDensity(GOdata.bac.enr, "GO:0006090", ranks=TRUE) # pyruvate metabolic process





# Enriched vs Unenriched
groups <- read.csv('groups.csv')
groups.rt <- groups[groups$RT_Enriched == 'True' | groups$RT_Unenriched == 'True', ]
geneList.rt.enr <- groups.rt$RT_Enriched == 'True' & groups.rt$RT_Unenriched == 'False'
geneList.rt.enr <- as.numeric(geneList.rt.enr)
names(geneList.rt.enr) <- groups.rt$X

sigFunc <- function(score) {return(score > 0.5)}

ontology <- 'MF'

GOdata.rt.enr <- new('topGOdata', ontology = ontology, allGenes = geneList.rt.enr, geneSel = sigFunc, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.rt.enr.fisher <- runTest(GOdata.rt.enr, statistic = 'fisher')
GenTable(GOdata.rt.enr, result.rt.enr.fisher, topNodes = 10)

ontology <- 'BP'

sigFunc <- function(score) {return(score > 0.5)}
GOdata.rt.enr <- new('topGOdata', ontology = ontology, allGenes = geneList.rt.enr, geneSel = sigFunc, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.rt.enr.fisher <- runTest(GOdata.rt.enr, statistic = 'fisher')
GenTable(GOdata.rt.enr, result.rt.enr.fisher, topNodes = 10)


# RAG Enr vs Unenr
groups.rag <- groups[groups$RAG_Enriched == 'True' | groups$RAG_Unenriched == 'True', ]
geneList.rag.enr <- groups.rag$RAG_Enriched == 'True' & groups.rag$RAG_Unenriched == 'False'
geneList.rag.enr <- as.numeric(geneList.rag.enr)
names(geneList.rag.enr) <- groups.rag$X


ontology <- 'MF'

sigFunc <- function(score) {return(score > 0.5)}
GOdata.rag.enr <- new('topGOdata', ontology = ontology, allGenes = geneList.rag.enr, geneSel = sigFunc, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.rag.enr.fisher <- runTest(GOdata.rag.enr, statistic = 'fisher')
GenTable(GOdata.rag.enr, result.rag.enr.fisher, topNodes = 10)

ontology <- 'BP'

sigFunc <- function(score) {return(score > 0.5)}
GOdata.rag.enr <- new('topGOdata', ontology = ontology, allGenes = geneList.rag.enr, geneSel = sigFunc, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.rag.enr.fisher <- runTest(GOdata.rag.enr, statistic = 'fisher')
GenTable(GOdata.rag.enr, result.rag.enr.fisher, topNodes = 10)



#Unenriched vs Enriched
groups <- read.csv('groups.csv')
geneList.rt.un <- groups.rt$RT_Enriched == 'False' & groups.rt$RT_Unenriched == 'True'
geneList.rt.un <- as.numeric(geneList.rt.un)
names(geneList.rt.un) <- groups.rt$X

ontology <- 'MF'

sigFunc <- function(score) {return(score > 0.5)}
GOdata.rt.un <- new('topGOdata', ontology = ontology, allGenes = geneList.rt.un, geneSel = sigFunc, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.rt.un.fisher <- runTest(GOdata.rt.un, statistic = 'fisher')
GenTable(GOdata.rt.un, result.rt.un.fisher, topNodes = 10)

ontology <- 'BP'

sigFunc <- function(score) {return(score > 0.5)}
GOdata.rt.un <- new('topGOdata', ontology = ontology, allGenes = geneList.rt.un, geneSel = sigFunc, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.rt.un.fisher <- runTest(GOdata.rt.un, statistic = 'fisher')
GenTable(GOdata.rt.un, result.rt.un.fisher, topNodes = 10)



# RAG Unenriched vs Enriched
geneList.rag.un <- groups.rag$RAG_Enriched == 'False' & groups.rag$RAG_Unenriched == 'True'
geneList.rag.un <- as.numeric(geneList.rag.un)
names(geneList.rag.un) <- groups.rag$X


ontology <- 'MF'

sigFunc <- function(score) {return(score > 0.5)}
GOdata.rag.un <- new('topGOdata', ontology = ontology, allGenes = geneList.rag.un, geneSel = sigFunc, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.rag.un.fisher <- runTest(GOdata.rag.un, statistic = 'fisher')
GenTable(GOdata.rag.un, result.rag.un.fisher, topNodes = 10)

ontology <- 'BP'

sigFunc <- function(score) {return(score > 0.5)}
GOdata.rag.un <- new('topGOdata', ontology = ontology, allGenes = geneList.rag.un, geneSel = sigFunc, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.rag.un.fisher <- runTest(GOdata.rag.un, statistic = 'fisher')
GenTable(GOdata.rag.un, result.rag.un.fisher, topNodes = 10)

