# 20190111_plot_kraken_res.r

# ------ PRELIMINARY SETTINGS ------

# import colours to remain consistent between all plots
source("/Users/Clara/Projects/colors/colors.R")
source("/Users/Clara/Projects/colors/colors2.R")

# set the directory to save img
img.path = "/Users/Clara/Projects/diary/graphics/plots/"
library(stringr)
today <- str_replace_all(Sys.Date(), "-", "") # concatenate the date following the format "yyyymmdd"
# the date is used to save output files such as plots without overwritting previous work (on a daily basis)

# ------ LIBRARIES ------

library(ggplot2)
library(tikzDevice)
library(ggpubr)
library(ghibli) # for pretty colours
library(phyloseq)

# ----- IMPORT DATA -----
# Input files are biom.json files that were generated from Kraken2 (with standard db)

path <- "/Users/Clara/Projects/osd/res/kraken2/miseq"

# check the list of files
list.files(path)

list.files <- list.files(path, pattern = ".biom.json$", recursive = TRUE)
list.files

biom.files <- list.files[grepl(".biom.json$", list.files)]
biom.paths <- file.path(path, biom.files)
file.names <- sapply(strsplit(sapply(strsplit(biom.files, "[.]"), `[`, 1), "/", 1), `[`, 2)

# ----- METADATA -----
# Generate the metadata table from file names

meta.data <- data.frame(
	"file.name" = file.names,
	"ena.id" = sapply(strsplit(file.names, "_"), `[`, 1),
	"sequencer" = sapply(strsplit(biom.paths, "[/]"), `[`, 8),
	"library" = sapply(strsplit(biom.paths, "[/]"), `[`, 9),
	"smp.cpn" = NA, #"OSD-Jun-2014",
	"stn.id" = NA,
	"stn.num" = NA,
	"smp.date" = NA, #as.Date("20140620", "%Y%m%d"),
	"lat" = NA, # DD
	"lon" = NA,
	"smp.depth" = NA,
	"temp" = NA, # celcius degrees
	"salt" = NA, # psu
	"primer" = NA, 
	"kraken2.db" = sapply(strsplit(file.names, "_"), `[`, 2)
	)
rownames(meta.data) <- file.names 

meta.data[meta.data$file.name=="ERR771024_custom",]$"stn.id" <- "FX0"
meta.data[meta.data$file.name=="ERR771024_custom",]$lat <- 64.208333
meta.data[meta.data$file.name=="ERR771024_custom",]$lon <- -22.015
meta.data[meta.data$file.name=="ERR771024_custom",]$smp.depth <- 0
meta.data[meta.data$file.name=="ERR771024_custom",]$temp <- 11
meta.data[meta.data$file.name=="ERR771024_custom",]$salt <- 31.2

meta.data[meta.data$file.name=="ERR771024_silva",]$"stn.id" <- "FX0"
meta.data[meta.data$file.name=="ERR771024_silva",]$lat <- 64.208333
meta.data[meta.data$file.name=="ERR771024_silva",]$lon <- -22.015
meta.data[meta.data$file.name=="ERR771024_silva",]$smp.depth <- 0
meta.data[meta.data$file.name=="ERR771024_silva",]$temp <- 11
meta.data[meta.data$file.name=="ERR771024_silva",]$salt <- 31.2

meta.data[meta.data$file.name=="ERR867776_custom",]$"stn.id" <- "FX0"
meta.data[meta.data$file.name=="ERR867776_custom",]$lat <- 64.208333
meta.data[meta.data$file.name=="ERR867776_custom",]$lon <- -22.015
meta.data[meta.data$file.name=="ERR867776_custom",]$smp.depth <- 0
meta.data[meta.data$file.name=="ERR867776_custom",]$temp <- 11
meta.data[meta.data$file.name=="ERR867776_custom",]$salt <- 31.2
meta.data[meta.data$file.name=="ERR867776_custom",]$primer <- "16S"

meta.data[meta.data$file.name=="ERR867776_silva",]$"stn.id" <- "FX0"
meta.data[meta.data$file.name=="ERR867776_silva",]$lat <- 64.208333
meta.data[meta.data$file.name=="ERR867776_silva",]$lon <- -22.015
meta.data[meta.data$file.name=="ERR867776_silva",]$smp.depth <- 0
meta.data[meta.data$file.name=="ERR867776_silva",]$temp <- 11
meta.data[meta.data$file.name=="ERR867776_silva",]$salt <- 31.2
meta.data[meta.data$file.name=="ERR867776_silva",]$primer <- "16S"

meta.data[meta.data$file.name=="ERR867777_custom",]$"stn.id" <- "FX0"
meta.data[meta.data$file.name=="ERR867777_custom",]$lat <- 64.208333
meta.data[meta.data$file.name=="ERR867777_custom",]$lon <- -22.015
meta.data[meta.data$file.name=="ERR867777_custom",]$smp.depth <- 0
meta.data[meta.data$file.name=="ERR867777_custom",]$temp <- 11
meta.data[meta.data$file.name=="ERR867777_custom",]$salt <- 31.2
meta.data[meta.data$file.name=="ERR867777_custom",]$primer <- "18S"

meta.data[meta.data$file.name=="ERR867777_silva",]$"stn.id" <- "FX0"
meta.data[meta.data$file.name=="ERR867777_silva",]$lat <- 64.208333
meta.data[meta.data$file.name=="ERR867777_silva",]$lon <- -22.015
meta.data[meta.data$file.name=="ERR867777_silva",]$smp.depth <- 0
meta.data[meta.data$file.name=="ERR867777_silva",]$temp <- 11
meta.data[meta.data$file.name=="ERR867777_silva",]$salt <- 31.2
meta.data[meta.data$file.name=="ERR867777_silva",]$primer <- "18S"

meta.data[meta.data$file.name=="ERR867778_custom",]$"stn.id" <- "FX0"
meta.data[meta.data$file.name=="ERR867778_custom",]$lat <- 64.208333
meta.data[meta.data$file.name=="ERR867778_custom",]$lon <- -22.015
meta.data[meta.data$file.name=="ERR867778_custom",]$smp.depth <- 20
meta.data[meta.data$file.name=="ERR867778_custom",]$temp <- 11
meta.data[meta.data$file.name=="ERR867778_custom",]$salt <- 31.2
meta.data[meta.data$file.name=="ERR867778_custom",]$primer <- "16S"

meta.data[meta.data$file.name=="ERR867778_silva",]$"stn.id" <- "FX0"
meta.data[meta.data$file.name=="ERR867778_silva",]$lat <- 64.208333
meta.data[meta.data$file.name=="ERR867778_silva",]$lon <- -22.015
meta.data[meta.data$file.name=="ERR867778_silva",]$smp.depth <- 20
meta.data[meta.data$file.name=="ERR867778_silva",]$temp <- 11
meta.data[meta.data$file.name=="ERR867778_silva",]$salt <- 31.2
meta.data[meta.data$file.name=="ERR867778_silva",]$primer <- "16S"

meta.data[meta.data$file.name=="ERR867779_custom",]$"stn.id" <- "FX0"
meta.data[meta.data$file.name=="ERR867779_custom",]$lat <- 64.208333
meta.data[meta.data$file.name=="ERR867779_custom",]$lon <- -22.015
meta.data[meta.data$file.name=="ERR867779_custom",]$smp.depth <- 20
meta.data[meta.data$file.name=="ERR867779_custom",]$temp <- 11
meta.data[meta.data$file.name=="ERR867779_custom",]$salt <- 31.2
meta.data[meta.data$file.name=="ERR867779_custom",]$primer <- "18S"

meta.data[meta.data$file.name=="ERR867779_silva",]$"stn.id" <- "FX0"
meta.data[meta.data$file.name=="ERR867779_silva",]$lat <- 64.208333
meta.data[meta.data$file.name=="ERR867779_silva",]$lon <- -22.015
meta.data[meta.data$file.name=="ERR867779_silva",]$smp.depth <- 0
meta.data[meta.data$file.name=="ERR867779_silva",]$temp <- 10.1
meta.data[meta.data$file.name=="ERR867779_silva",]$salt <- 48
meta.data[meta.data$file.name=="ERR867779_silva",]$primer <- "18S"

# Eyjafjordur EFJ5
meta.data[meta.data$file.name=="ERR867894_custom",]$"stn.id" <- "EFJ5"
meta.data[meta.data$file.name=="ERR867894_custom",]$lat <- 66.1316
meta.data[meta.data$file.name=="ERR867894_custom",]$lon <- -18.7902 
meta.data[meta.data$file.name=="ERR867894_custom",]$smp.depth <- 20
meta.data[meta.data$file.name=="ERR867894_custom",]$temp <- 11
meta.data[meta.data$file.name=="ERR867894_custom",]$salt <- 31.2
meta.data[meta.data$file.name=="ERR867894_custom",]$primer <- "16S"

meta.data[meta.data$file.name=="ERR867894_silva",]$"stn.id" <- "EFJ5"
meta.data[meta.data$file.name=="ERR867894_silva",]$lat <- 66.1316
meta.data[meta.data$file.name=="ERR867894_silva",]$lon <- -18.7902 
meta.data[meta.data$file.name=="ERR867894_silva",]$smp.depth <- 0
meta.data[meta.data$file.name=="ERR867894_silva",]$temp <- 10.1
meta.data[meta.data$file.name=="ERR867894_silva",]$salt <- 48
meta.data[meta.data$file.name=="ERR867894_silva",]$primer <- "16S"

meta.data[meta.data$file.name=="ERR867936_custom",]$"stn.id" <- "IFJ0" # Isafjordur
meta.data[meta.data$file.name=="ERR867936_custom",]$lat <- 65.9449
meta.data[meta.data$file.name=="ERR867936_custom",]$lon <- -22.4192 
meta.data[meta.data$file.name=="ERR867936_custom",]$smp.depth <- 0
meta.data[meta.data$file.name=="ERR867936_custom",]$temp <- 7.6
meta.data[meta.data$file.name=="ERR867936_custom",]$salt <- 19
meta.data[meta.data$file.name=="ERR867936_custom",]$primer <- "16S"

meta.data[meta.data$file.name=="ERR867936_silva",]$"stn.id" <- "IFJ0"
meta.data[meta.data$file.name=="ERR867936_silva",]$lat <- 65.9449
meta.data[meta.data$file.name=="ERR867936_silva",]$lon <- -22.4192 
meta.data[meta.data$file.name=="ERR867936_silva",]$smp.depth <- 0
meta.data[meta.data$file.name=="ERR867936_silva",]$temp <- 7.6
meta.data[meta.data$file.name=="ERR867936_silva",]$salt <- 19
meta.data[meta.data$file.name=="ERR867936_silva",]$primer <- "16S"

meta.data[meta.data$file.name=="ERR867937_custom",]$"stn.id" <- "IFJ0" # Isafjordur
meta.data[meta.data$file.name=="ERR867937_custom",]$lat <- 65.9449
meta.data[meta.data$file.name=="ERR867937_custom",]$lon <- -22.4192 
meta.data[meta.data$file.name=="ERR867937_custom",]$smp.depth <- 0
meta.data[meta.data$file.name=="ERR867937_custom",]$temp <- 7.6
meta.data[meta.data$file.name=="ERR867937_custom",]$salt <- 19
meta.data[meta.data$file.name=="ERR867937_custom",]$primer <- "18S"

meta.data[meta.data$file.name=="ERR867937_silva",]$"stn.id" <- "IFJ0"
meta.data[meta.data$file.name=="ERR867937_silva",]$lat <- 65.9449
meta.data[meta.data$file.name=="ERR867937_silva",]$lon <- -22.4192 
meta.data[meta.data$file.name=="ERR867937_silva",]$smp.depth <- 0
meta.data[meta.data$file.name=="ERR867937_silva",]$temp <- 7.6
meta.data[meta.data$file.name=="ERR867937_silva",]$salt <- 19
meta.data[meta.data$file.name=="ERR867937_silva",]$primer <- "18S"

meta.data[meta.data$file.name=="ERR867938_custom",]$"stn.id" <- "IFJ0" # Isafjordur
meta.data[meta.data$file.name=="ERR867938_custom",]$lat <- 65.9449
meta.data[meta.data$file.name=="ERR867938_custom",]$lon <- -22.4192 
meta.data[meta.data$file.name=="ERR867938_custom",]$smp.depth <- 15
meta.data[meta.data$file.name=="ERR867938_custom",]$temp <- 7.5
meta.data[meta.data$file.name=="ERR867938_custom",]$salt <- 29.2
meta.data[meta.data$file.name=="ERR867938_custom",]$primer <- "18S"

meta.data[meta.data$file.name=="ERR867938_silva",]$"stn.id" <- "IFJ0"
meta.data[meta.data$file.name=="ERR867938_silva",]$lat <- 65.9449
meta.data[meta.data$file.name=="ERR867938_silva",]$lon <- -22.4192 
meta.data[meta.data$file.name=="ERR867938_silva",]$smp.depth <- 15
meta.data[meta.data$file.name=="ERR867938_silva",]$temp <- 7.5
meta.data[meta.data$file.name=="ERR867938_silva",]$salt <- 29.2
meta.data[meta.data$file.name=="ERR867938_silva",]$primer <- "18S"

meta.data[meta.data$file.name=="ERR867939_custom",]$"stn.id" <- "IFJ0"
meta.data[meta.data$file.name=="ERR867939_custom",]$lat <- 65.9449
meta.data[meta.data$file.name=="ERR867939_custom",]$lon <- -22.4192
meta.data[meta.data$file.name=="ERR867939_custom",]$smp.depth <- 15
meta.data[meta.data$file.name=="ERR867939_custom",]$temp <- 7.5
meta.data[meta.data$file.name=="ERR867939_custom",]$salt <- 29.2
meta.data[meta.data$file.name=="ERR867939_custom",]$primer <- "16S"

meta.data[meta.data$file.name=="ERR867939_silva",]$"stn.id" <- "IFJ0"
meta.data[meta.data$file.name=="ERR867939_silva",]$lat <- 65.9449
meta.data[meta.data$file.name=="ERR867939_silva",]$lon <- -22.4192
meta.data[meta.data$file.name=="ERR867939_silva",]$smp.depth <- 15
meta.data[meta.data$file.name=="ERR867939_silva",]$temp <- 7.5
meta.data[meta.data$file.name=="ERR867939_silva",]$salt <- 29.2
meta.data[meta.data$file.name=="ERR867939_silva",]$primer <- "16S"

meta.data[meta.data$file.name=="ERR771077_custom",]$"stn.id" <- "EFJ5"
meta.data[meta.data$file.name=="ERR771077_custom",]$lat <- 66.1316 
meta.data[meta.data$file.name=="ERR771077_custom",]$lon <- -18.7902
meta.data[meta.data$file.name=="ERR771077_custom",]$smp.depth <- 0
meta.data[meta.data$file.name=="ERR771077_custom",]$temp <- 10.1
meta.data[meta.data$file.name=="ERR771077_custom",]$salt <- 48
meta.data[meta.data$file.name=="ERR771077_custom",]$primer <- NA

meta.data[meta.data$file.name=="ERR771077_silva",]$"stn.id" <- "EFJ5"
meta.data[meta.data$file.name=="ERR771077_silva",]$lat <- 66.1316 
meta.data[meta.data$file.name=="ERR771077_silva",]$lon <- -18.7902
meta.data[meta.data$file.name=="ERR771077_silva",]$smp.depth <- 0
meta.data[meta.data$file.name=="ERR771077_silva",]$temp <- 10.1
meta.data[meta.data$file.name=="ERR771077_silva",]$salt <- 48
meta.data[meta.data$file.name=="ERR771077_silva",]$primer <- NA

meta.data[meta.data$file.name=="ERR771098_custom",]$"stn.id" <- "IFJ0"
meta.data[meta.data$file.name=="ERR771098_custom",]$lat <- 65.9449
meta.data[meta.data$file.name=="ERR771098_custom",]$lon <- -22.4192
meta.data[meta.data$file.name=="ERR771098_custom",]$smp.depth <- 15
meta.data[meta.data$file.name=="ERR771098_custom",]$temp <- 7.5
meta.data[meta.data$file.name=="ERR771098_custom",]$salt <- 29.2
meta.data[meta.data$file.name=="ERR771098_custom",]$primer <- NA

meta.data[meta.data$file.name=="ERR771098_silva",]$"stn.id" <- "IFJ0"
meta.data[meta.data$file.name=="ERR771098_silva",]$lat <- 65.9449
meta.data[meta.data$file.name=="ERR771098_silva",]$lon <- -22.4192
meta.data[meta.data$file.name=="ERR771098_silva",]$smp.depth <- 15
meta.data[meta.data$file.name=="ERR771098_silva",]$temp <- 7.5
meta.data[meta.data$file.name=="ERR771098_silva",]$salt <- NA
meta.data[meta.data$file.name=="ERR771098_silva",]$primer <- NA


# meta.data[meta.data$file.name=="ERR867938_custom",]$"stn.id" <- NA
# meta.data[meta.data$file.name=="ERR867938_custom",]$lat <- NA
# meta.data[meta.data$file.name=="ERR867938_custom",]$lon <- NA
# meta.data[meta.data$file.name=="ERR867938_custom",]$smp.depth <- NA
# meta.data[meta.data$file.name=="ERR867938_custom",]$temp <- NA
# meta.data[meta.data$file.name=="ERR867938_custom",]$salt <- NA
# meta.data[meta.data$file.name=="ERR867938_custom",]$primer <- NA

# meta.data[meta.data$file.name=="ERR867938_silva",]$"stn.id" <- NA
# meta.data[meta.data$file.name=="ERR867938_silva",]$lat <- NA
# meta.data[meta.data$file.name=="ERR867938_silva",]$lon <- NA
# meta.data[meta.data$file.name=="ERR867938_silva",]$smp.depth <- NA
# meta.data[meta.data$file.name=="ERR867938_silva",]$temp <- NA
# meta.data[meta.data$file.name=="ERR867938_silva",]$salt <- NA
# meta.data[meta.data$file.name=="ERR867938_silva",]$primer <- NA


# ------ DATASETS ------

bioms <- list()
for(i in seq_along(biom.paths)) {
	bioms[i] <- import_biom(biom.paths[i], 
 	treefilename=NULL, refseqfilename=NULL, refseqFunction=readDNAStringSet, refseqArgs=NULL,
 	parseFunction=parse_taxonomy_greengenes, version=1.0)
 	colnames(otu_table(bioms[[i]])) <- file.names[i]
}

# data <- merge_phyloseq(bioms[[1]], bioms[[2]], bioms[[3]], bioms[[4]], bioms[[5]], bioms[[6]], bioms[[7]], bioms[[8]], bioms[[9]])

# Create a big phyloseq object with all json biom 
data <- bioms[[1]]
for (i in seq(2, length(bioms))) {
	data <- merge_phyloseq(data, bioms[[i]])
}
sample_data(data) <- meta.data

# ----- REMOVE WEIRD DATA ----

data <- subset_taxa(data, Kingdom!="Holozoa")
data <- subset_taxa(data, Kingdom!="Nucletmycea")

# ----- VISUALS ------

# collapse taxonomy to Kingdom rank
data.kingdoms <- tax_glom(data, taxrank=rank_names(data)[1], 
         NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

p = plot_bar(data.kingdoms, x="library", fill = "Kingdom") +
	geom_bar( aes(color = Kingdom, fill = Kingdom ), 
 		stat = 'identity', 
 		position = 'stack') 

p + facet_wrap(~kraken2.db)

# ------ BACTERIA ----- 
bakt <- subset_taxa(data, Kingdom=="Bacteria")

# collapse taxonomy to phylum rank
p.bakt <- tax_glom(bakt, 
  taxrank=rank_names(bakt)[2], 
  NArm=TRUE, 
  bad_empty=c(NA, "", " ", "\t"))

p = plot_bar(p.bakt, x="library", fill = "Phylum") +
	geom_bar( aes(color = Phylum, fill = Phylum ), 
 		stat = 'identity', 
 		position = 'stack') 

p + facet_wrap(~kraken2.db)

c.data <- subset_samples(p.bakt, kraken2.db=="silva")
rownames(sample_data(c.data)) <- 

rownames(meta.data) <- file.names 

p = plot_bar(c.data, x="library", fill = "Phylum") +
	geom_bar( aes(color = Phylum, fill = Phylum ), 
 		stat = 'identity', 
 		position = 'stack') 

p 

# ------ ALPHA DIVERSITY ------
# In ecology, alpha diversity (α-diversity) 
# is the mean species diversity in sites or 
# habitats at a local scale.

# http://www.metagenomics.wiki/pdf/definition/alpha-beta-diversity

pAlpha = plot_richness(c.data,
   shape = "library",
   color = "stn.id",
   measures = c("Observed", "Chao1", "Shannon", "InvSimpson"),
   title = "OSD Alpha Diversity")
pAlpha + geom_point(size = 5)

# here it is meaningless to use beta-diversity
# because we compare the same sample!

# ----- RAREFACTION -----

r.data <- rarefy_even_depth(data, replace=TRUE)

kr.data <- tax_glom(r.data, taxrank=rank_names(r.data)[1], 
         NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

p = plot_bar(kr.data, x="library", fill = "Kingdom") +
	geom_bar( aes(color = Kingdom, fill = Kingdom ), 
 		stat = 'identity', 
 		position = 'stack') 

p.alpha = plot_richness(kr.data,
   shape = "library",
   color = "stn.id",
   measures = c("Observed", "Chao1", "Shannon", "InvSimpson"),
   title = "OSD Alpha Diversity")
p.alpha + geom_point(size = 5)


bakt <- subset_taxa(data, Kingdom=="Bacteria")

r.bakt <- rarefy_even_depth(bakt, replace=TRUE)

p = plot_bar(r.bakt, x="library", fill = "Phylum") +
	geom_bar( aes(color = Phylum, fill = Phylum ), 
 		stat = 'identity', 
 		position = 'stack') 

p.alpha = plot_richness(r.bakt,
   shape = "library",
   color = "stn.id",
   measures = c("Observed", "Chao1", "Shannon", "InvSimpson"),
   title = "OSD Alpha Diversity")
p.alpha + geom_point(size = 5)

# ---- DESeq -----

library("DESeq2")
packageVersion("DESeq2")

kostic <- prune_samples(sample_sums(data) > 500, data)
kostic

diagdds = phyloseq_to_deseq2(kostic, ~ library)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(kostic)[rownames(sigtab), ], "matrix"))
head(sigtab)


posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
#posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

# Plot Results

library("ggplot2")
theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Genus))
# Kingdom order
x = tapply(sigtabgen$log2FoldChange, 
	sigtabgen$Kingdom, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Kingdom = factor(as.character(sigtabgen$Kingdom), 
	levels=names(x))

# Phylum order
x = tapply(sigtabgen$log2FoldChange, 
	sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
ggplot(sigtabgen, aes(y=Phylum, x=log2FoldChange, color=Kingdom)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))



