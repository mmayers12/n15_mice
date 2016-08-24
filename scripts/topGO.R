library(topGO)
library(stats)
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


sampleTable <- read.csv(file.path(data.dir, 'filt_metadata.csv'), row.names=1)
sampleTable$technical <- as.logical(sampleTable$technical)
locusTable <- read.csv(file.path(data.dir, 'loci_annot.csv'))


# Start by looking at human and mouse proteins only
mh.prots <- locusTable[locusTable$mouse_human == 'True',]$X
mh.pvals <- unenr_pvals[unenr_pvals$X %in% mh.prots, ]
geneList.mh <- mh.pvals$p_value
geneList.mh <- mh.pvals$adjFC
names(geneList.mh) <- mh.pvals$X



mh.pvals.up <- mh.pvals[mh.pvals$log2FoldChange >= 0, ]
geneList.mh.up <- mh.pvals.up$p_value
names(geneList.mh.up) <- mh.pvals.up$X

mh.pvals.down <- mh.pvals[mh.pvals$log2FoldChange < 0, ]
geneList.mh.down <- mh.pvals.down$p_value
names(geneList.mh.down) <- mh.pvals.down$X

geneList.mh.updown <- c(geneList.mh.down * -1, geneList.mh.up)


ontology <- 'MF'
topDiffGenes <- function(score) {return(score < 0.01)}

geneList.mh  <- -1*geneList.mh

GOdata.mh <- new('topGOdata', ontology = ontology, allGenes = geneList.mh, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.mh.fisher <- runTest(GOdata.mh, statistic = 'fisher')
result.mh.ks <- runTest(GOdata.mh, statistic = 'ks')
GenTable(GOdata.mh, result.mh.fisher, topNodes = 10)
GenTable(GOdata.mh, result.mh.ks, topNodes = 10)

GOdata.mh.up <- new('topGOdata', ontology = ontology, allGenes = geneList.mh.up, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.mh.up.fisher <- runTest(GOdata.mh.up, statistic = 'fisher')
result.mh.up.ks <- runTest(GOdata.mh.up, statistic = 'ks')
GenTable(GOdata.mh.up, result.mh.up.fisher, topNodes = 10)
GenTable(GOdata.mh.up, result.mh.up.ks, topNodes = 10)

GOdata.mh.down <- new('topGOdata', ontology = ontology, allGenes = geneList.mh.down, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.mh.down.fisher <- runTest(GOdata.mh.down, statistic = 'fisher')
result.mh.down.ks <- runTest(GOdata.mh.down, statistic = 'ks')
GenTable(GOdata.mh.down, result.mh.down.fisher, topNodes = 10)
GenTable(GOdata.mh.down, result.mh.down.ks, topNodes = 10)


geneList.mh.updown <- -1*geneList.mh.updown

GOdata.mh.updown <- new('topGOdata', ontology = ontology, allGenes = geneList.mh.updown, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.mh.updown.fisher <- runTest(GOdata.mh.updown, statistic = 'fisher')
result.mh.updown.ks <- runTest(GOdata.mh.updown, statistic = 'ks')
GenTable(GOdata.mh.updown, result.mh.updown.fisher, topNodes = 10)
GenTable(GOdata.mh.updown, result.mh.updown.ks, topNodes = 10)





# Now for bacterial proteins only
tmp <- locusTable$mouse_human == 'False' & locusTable$lca != '35823' # Filtering out anthrospira differences 35823
tmp[is.na(tmp)] <- TRUE
bac.prots <- locusTable[tmp, ]$X
bac.pvals <- unenr_pvals[unenr_pvals$X %in% bac.prots, ]
geneList.bac <- bac.pvals$padj
geneList.bac <- bac.pvals$adjFC
names(geneList.bac) <- bac.pvals$X

bac.pvals.up <- bac.pvals[bac.pvals$log2FoldChange > 0, ]
geneList.bac.up <- bac.pvals.up$p_value
names(geneList.bac.up) <- bac.pvals.up$X

bac.pvals.down <- bac.pvals[bac.pvals$log2FoldChange < 0, ]
geneList.bac.down <- bac.pvals.down$p_value
names(geneList.bac.down) <- bac.pvals.down$X

geneList.bac.updown <- c(geneList.bac.down * -1, geneList.bac.up)


## Look at what's inside 
ontology <- 'MF'

geneList.bac <- -1*geneList.bac

GOdata.bac <- new('topGOdata', ontology = ontology, allGenes = geneList.bac, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.bac.fisher <- runTest(GOdata.bac, statistic = 'fisher')
result.bac.ks <- runTest(GOdata.bac, statistic = 'ks')
GenTable(GOdata.bac, result.bac.fisher, topNodes = 10)
GenTable(GOdata.bac, result.bac.ks, topNodes = 10)

GOdata.bac.up <- new('topGOdata', ontology = ontology, allGenes = geneList.bac.up, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.bac.up.fisher <- runTest(GOdata.bac.up, statistic = 'fisher')
result.bac.up.ks <- runTest(GOdata.bac.up, statistic = 'ks')
GenTable(GOdata.bac.up, result.bac.up.fisher, topNodes = 10)
GenTable(GOdata.bac.up, result.bac.up.ks, topNodes = 10)

GOdata.bac.down <- new('topGOdata', ontology = ontology, allGenes = geneList.bac.down, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.bac.down.fisher <- runTest(GOdata.bac.down, statistic = 'fisher')
result.bac.down.ks <- runTest(GOdata.bac.down, statistic = 'ks')
GenTable(GOdata.bac.down, result.bac.down.fisher, topNodes = 10)
GenTable(GOdata.bac.down, result.bac.down.ks, topNodes = 10)

geneList.bac.updown <- -1 * geneList.bac.updown

GOdata.bac.updown <- new('topGOdata', ontology = ontology, allGenes = geneList.bac.updown, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.bac.updown.fisher <- runTest(GOdata.bac.updown, statistic = 'fisher')
result.bac.updown.ks <- runTest(GOdata.bac.updown, statistic = 'ks')
GenTable(GOdata.bac.updown, result.bac.updown.ks, topNodes = 10)


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
ontology <- 'MF'

# Tcell vs RAG
geneList.bac.enr <- -1*geneList.bac.enr # KS goes from - to +
                                        # do this to put upreg in Tcell first

GOdata.bac <- new('topGOdata', ontology = ontology, allGenes = geneList.bac, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.bac.ks <- runTest(GOdata.bac, statistic = 'ks')
GenTable(GOdata.bac, result.bac.ks, topNodes = 10)

ontology <- 'BP'
GOdata.bac <- new('topGOdata', ontology = ontology, allGenes = geneList.bac, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.bac.ks <- runTest(GOdata.bac, statistic = 'ks')
GenTable(GOdata.bac, result.bac.ks, topNodes = 10)

# RAG vs TCell
ontology <- 'MF'
geneList.bac.enr <- -1*geneList.bac.enr 

GOdata.bac <- new('topGOdata', ontology = ontology, allGenes = geneList.bac, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.bac.ks <- runTest(GOdata.bac, statistic = 'ks')
GenTable(GOdata.bac, result.bac.ks, topNodes = 10)

ontology <- 'BP'
GOdata.bac <- new('topGOdata', ontology = ontology, allGenes = geneList.bac, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = gene2GOID)
result.bac.ks <- runTest(GOdata.bac, statistic = 'ks')
GenTable(GOdata.bac, result.bac.ks, topNodes = 10)






# Enriched vs Unenriched
groups <- read.csv('groups.csv')
groups.rt <- groups[groups$RT_Enriched == 'True' | groups$RT_Unenriched == 'True', ]
geneList.rt.enr <- groups.rt$RT_Enriched == 'True' & groups.rt$RT_Unenriched == 'False'
geneList.rt.enr <- as.numeric(geneList.rt.enr)
names(geneList.rt.enr) <- groups.rt$X

ontology <- 'MF'

sigFunc <- function(score) {return(score > 0.5)}
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

