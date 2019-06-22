library('biomaRt')
library('xml2')
library("illuminaHumanv4.db")
library('lumi')
library(lumiHumanIDMapping)

listMarts()
ensembl=useMart("ensembl")

datasets <- listDatasets(ensembl)
head(datasets)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

# filters
filters = listFilters(ensembl)
filters[1:5,]

# attributes
attributes = listAttributes(ensembl)
attributes[1:5,]


###### illumina run ######
# searchAttributes(mart = ensembl, pattern = "illum")
# searchAttributes(mart = ensembl, pattern = "ensembl.*id")

# searchFilters(mart = ensembl, pattern = "illum")
# searchFilters(mart = ensembl, pattern = "array")

getwd()
setwd('~/Documents/Masters/RNA_seq_classifier/Data')
illumina <- read.table('ill_probe.csv', sep = ',', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)
head(illumina)
length(unique(illumina$Probe_Id))
length(unique(illumina$Array_Address_Id))


setwd('/home/patrick/Code/R')
load('esets.RData')
setwd('~/Documents/Masters/RNA_seq_classifier/Data/Ciber_sort')

idx <- status['most_general'] == 'bacterial' |
  status['most_general'] == 'viral' |
  status['most_general'] == 'greyb' |
  status['most_general'] == 'greyv'|
  status['most_general'] == 'greyu'
sum(idx)

e.set.f <- as.data.frame(e.set[,idx])
dim(e.set.f)

# rownames(e.set.f)[1:10]

e.set.f$Probe_Id <- illumina$Probe_Id[match(rownames(e.set.f),illumina$Array_Address_Id)]
e.set.f[1:5,(ncol(e.set.f)-4):ncol(e.set.f)]


# testing
IlluminaID2nuID(e.set.f[1:5,(ncol(e.set.f)-4):ncol(e.set.f)][5], lib.mapping='lumiHumanIDMapping', chipVersion = NULL)

# Search_Key    ILMN_Gene Accession     Symbol Probe_Id       Array_Address_Id nuID                
# 6450255 "NM_182762.2" "7A5"     "NM_182762.2" "7A5"  "ILMN_1762337" "0006450255"     "Ku8QhfS0n_hIOABXuE"
# 2570615 "NM_130786.2" "A1BG"    "NM_130786.2" "A1BG" "ILMN_2055271" "0002570615"     "fqPEquJRRlSVSfL.8A"
# 6370619 "NM_130786.2" "A1BG"    "NM_130786.2" "A1BG" "ILMN_1736007" "0006370619"     "ckiehnugOno9d7vf1Q"
# 2600039 "NM_138932.1" "A1CF"    "NM_138932.1" "A1CF" "ILMN_2383229" "0002600039"     "x57Vw5B5Fbt5JUnQkI"
# 2650615 "NM_138933.1" "A1CF"    "NM_014576.2" "A1CF" "ILMN_1806310" "0002650615"     "ritxUH.kuHlYqjozpE"

nuID2IlluminaID('Ku8QhfS0n_hIOABXuE')
nuID2probeID('Ku8QhfS0n_hIOABXuE')

nuID2EntrezID('Ku8QhfS0n_hIOABXuE', lib.mapping='lumiHumanIDMapping')
nuID2RefSeqID('Ku8QhfS0n_hIOABXuE', lib.mapping='lumiHumanIDMapping')
nuID2targetID('Ku8QhfS0n_hIOABXuE', lib.mapping='lumiHumanIDMapping')



# real
con <- IlluminaID2nuID(e.set.f$Probe_Id, lib.mapping='lumiHumanIDMapping', chipVersion = NULL)
con.df <- as.data.frame(con)
dim(con.df)
head(con.df)

# check the e.set.f probe id match that in the conversion df
sum(con.df$Probe_Id == e.set.f$Probe_Id)

# check the conversion df if gene and symbol are equivalent = FALSE
sum(as.character(con.df$ILMN_Gene) == as.character(con.df$Symbol))

# reads in the mixture2 file used by Pre to check against
mixture2 <- read.table('mixture2.txt', sep = '\t', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)
dim(mixture2)

sum(mixture2$Gene == as.character(con.df$ILMN_Gene))
sum(mixture2$Gene == as.character(con.df$Symbol))

length(intersect(mixture2$Gene, as.character(con.df$ILMN_Gene)))
length(intersect(mixture2$Gene, as.character(con.df$Symbol)))
# looks like the gene name is closer to that used by pre in the mixture2

which(mixture2$Gene != as.character(con.df$ILMN_Gene))
mixture2$Gene[which(mixture2$Gene != as.character(con.df$ILMN_Gene))]
con.df$ILMN_Gene[which(mixture2$Gene != as.character(con.df$ILMN_Gene))]

ill.probes <- e.set.f$Probe_Id[which(mixture2$Gene != as.character(con.df$ILMN_Gene))]
ill.probes
length(ill.probes)
# ill.probes <- e.set.f$illumina_probe[1:15000]

df.1 <- getBM(attributes=c('illumina_humanht_12_v4', 'hgnc_symbol', 'ensembl_gene_id'), 
              filters = 'illumina_humanht_12_v4', 
              values = ill.probes, 
              mart = ensembl)
View(df.1)


# looks like when we do the biomaRt conversion that the hgnc gene names are capitalized just like our con.df ilmn_gene
# however we have some emptyy rows in our biomaRt conversion
# double checking these in our con.df conversion we can see they are not empty

df.1$illumina_humanht_12_v4[df.1$hgnc_symbol == '']
match(df.1$illumina_humanht_12_v4[df.1$hgnc_symbol == ''], con.df$Probe_Id)
con.df[match(df.1$illumina_humanht_12_v4[df.1$hgnc_symbol == ''], con.df$Probe_Id),]

dim(e.set.f)
dim(con.df)
e.set.f$Gene <- con.df$ILMN_Gene
View(cbind(Gene=con.df$ILMN_Gene, e.set.f))
View(cbind(Gene=mixture2$Gene, e.set.f))

# rip the gene column from the mixture2 and stick it on my version as col1
e.set.f <- cbind(Gene=mixture2$Gene, e.set.f)
dim(e.set.f)
View(head(e.set.f))

# saving and reading 
# setwd('/home/patrick/Documents/Masters/RNA_seq_classifier/Data/Ciber_sort')
# write.table(mixture2, file = "mixture2_rep1.txt", sep='\t', row.names=FALSE) # FAILED RUN - no overlapping genes
# write.table(mixture2, file = "mixture2_rep2.txt", sep='\t', row.names=TRUE) # 
dim(mixture2[,idx])
mixture2[,idx][1:5,1:5]


write.table(mixture2, file = "mixture2_rep2.txt", sep='\t', row.names=FALSE, quote = FALSE) # worked!
write.table(e.set.f, file = "eset1.txt", sep='\t', row.names=FALSE, quote = FALSE) # worked!!!!!!!!!!

















# # check that mixture2 genes are contained in e.set.f genes
# intersect(mixture2$Gene, e.set.f$Gene)
# 
# colnames(e.set.f)
# col.names <- c('Gene', colnames(e.set.f))
# # e.names <-c('Gene', cols.names)
# cols.reorder <- col.names[-length(col.names)]
# 
mix <- e.set.f
# # mix <- e.set.f[1:1000, cols.reorder]
# dim(mix)
# rownames(mix) <- NULL
# View(mix)
# 
# intersect(mixture2$Gene, mix$Gene)




mix_one <- read.table('mix_one.txt', sep = '\t', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)
LM22 <- read.table('LM22.txt', sep = '\t', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)
examp.mix <- read.table('ExampleMixtures-GEPs.txt', sep = '\t', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)
mixture2 <- read.table('mixture2.txt', sep = '\t', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)


dim(examp.mix)
dim(mixture2)
dim(mix_one)
dim(LM22)

LM22$Gene.symbol == examp.mix$GeneSymbol
LM22$Gene.symbol
examp.mix$GeneSymbol

intersect(LM22$Gene.symbol, examp.mix$GeneSymbol)
intersect(LM22$Gene.symbol, mixture2$Gene)
intersect(LM22$Gene.symbol, mix_one$Gene)

# its the same fucking list!
View(cbind(mix_one$Gene, mixture2$Gene))

mix_one$Gene <- mixture2$Gene

# setwd('/home/patrick/Documents/Masters/RNA_seq_classifier/Data/Ciber_sort')
# write.table(mix, file = "mix_one.txt", sep='\t', row.names=FALSE)



# rownames(mix_test)
# View(mix_test)



intersect(examp.mix$GeneSymbol, mixture2$Gene)

# pre mixture
mix2 <- read.table('mixture2.txt', sep = '\t', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)
dim(mix2)
rownames(mix2)
View(mix2[1:100,])

LM22.ref <- read.table('LM22-ref-sample.txt', sep = '\t', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)
dim(LM22.ref)
View(LM22.ref)


###############################################################################################
# ill.probs <- illumina$Probe_Id
# ill.probs <- illumina$Probe_Id[1:20000]
# ill.probs <- illumina$Probe_Id[10000:19999]
# ill.probs <- illumina$Probe_Id[40000:47323]
# ill.probs
# length(ill.probs)
# 
# df.1 <- getBM(attributes=c('illumina_humanht_12_v4', 'hgnc_symbol'), 
#       filters = 'illumina_humanht_12_v4', 
#       values = ill.probs, 
#       mart = ensembl)
# 
# df.1
# View(df.1)
# dim(df.1)
# length(unique(df.1$illumina_humanht_12_v4))
# length(unique(df.1$hgnc_symbol))



# ill_hugo_10000 <- read.table('ill_hugo_10000.csv', sep = ',', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)
# dim(ill_hugo_10000)



mysql  --user=genome --host=genome-mysql.cse.ucsc.edu -A -D hg19
select distinct G.gene,N.value from ensGtp as G, ensemblToGeneName as N where G.transcript=N.name and G.gene in ("ENSG00000183742", "ENSG00000121410") ;












# 
# LM22_ref <- read.table('LM22-ref-sample.txt', sep = '\t', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)
# View(LM22_ref)
# lm<-LM22_ref$Relabel
# length(lm)
# 
# df.3 <- getBM(attributes=c('hgnc_symbol', 'illumina_humanht_12_v4'), 
#               filters = 'hgnc_symbol', 
#               values = lm, 
#               mart = ensembl)
# dim(df.3)
# 
# View(df.3)
# 


