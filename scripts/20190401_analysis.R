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

ena.osd <- read.csv("/Users/Clara/Projects/osd/data/osd.csv")

meta.data <- data.frame(
  "file.name" = file.names,
  "ena.id" = sapply(strsplit(file.names, "_"), `[`, 1),
  "sequencer" = sapply(strsplit(biom.paths, "[/]"), `[`, 8),
  "library" = sapply(strsplit(biom.paths, "[/]"), `[`, 9),
  "smp.cpn" = "OSD-Jun-2014",
  "stn.id" = NA,
  "lat" = NA,
  "lon" = NA,
  "smp.depth" = NA,
  "temp" = NA,
  "salt" = NA,
  "primer" = NA,
  
  "kraken2.db" = sapply(strsplit(file.names, "_"), `[`, 2)
)

for (i in seq(1,dim(meta.data)[1],1)) {
  #print(meta.data[i,]$ena.id)
  for (j in seq(1,dim(ena.osd)[1],1)) {
    #print(ena.osd[j,]$ena.id)
    if(meta.data[i,]$ena.id == ena.osd[j,]$ena.id){
      meta.data[i,]$stn.id <- as.character(ena.osd[j,]$stn.id)
      meta.data[i,]$lat <- as.numeric(ena.osd[j,]$lat)
      meta.data[i,]$lon <- as.numeric(ena.osd[j,]$lon)
      meta.data[i,]$smp.depth <- as.numeric(ena.osd[j,]$smp.depth)
      meta.data[i,]$temp <- as.numeric(ena.osd[j,]$temp)
      meta.data[i,]$salt <- as.numeric(ena.osd[j,]$salt)
      meta.data[i,]$primer <- as.character(ena.osd[j,]$primer)
    }
  }
}

rownames(meta.data) <- file.names

# ------ DATASETS ------

bioms <- list()
for(i in seq_along(biom.paths)) {
  bioms[i] <- import_biom(biom.paths[i], 
                          treefilename=NULL, refseqfilename=NULL, refseqFunction=readDNAStringSet, refseqArgs=NULL,
                          parseFunction=parse_taxonomy_greengenes, version=1.0)
  colnames(otu_table(bioms[[i]])) <- file.names[i]
}

# Create a big phyloseq object with all json biom 
data <- bioms[[1]]
for (i in seq(2, length(bioms))) {
  data <- merge_phyloseq(data, bioms[[i]])
}
sample_data(data) <- meta.data

# ----- REMOVE WEIRD DATA ----

data <- subset_taxa(data, Kingdom!="Holozoa")
data <- subset_taxa(data, Kingdom!="Nucletmycea")
data <- subset_taxa(data, Genus!="Homo")

# collapse taxonomy to Kingdom rank
data.kingdoms <- tax_glom(data, taxrank=rank_names(data)[1], 
                          NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

p = plot_bar(data.kingdoms, x="library", fill = "Kingdom") +
  geom_bar( aes(color = Kingdom, fill = Kingdom ), 
            stat = 'identity', 
            position = 'stack') +
  theme_pubclean()

pdf(file = paste(img.path, today, "_all_data_db_comparison.pdf", sep="", collapse=NULL),  width=6, height=6)
p + facet_wrap(~kraken2.db)
dev.off()

# ----- ONE SAMPLE -----

fx0 <- subset_samples(data, stn.id=="FX0" & smp.depth == 0)

# ----- COMPARISON DB -----

# collapse taxonomy to Kingdom rank
data.kingdoms <- tax_glom(fx0, taxrank=rank_names(fx0)[1], 
                          NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

p = plot_bar(data.kingdoms, x="library", fill = "Kingdom") +
  geom_bar( aes(color = Kingdom, fill = Kingdom ), 
            stat = 'identity', 
            position = 'stack') +
  theme_pubclean()

pdf(file = paste(img.path, today, "_fx0_db_comparison.pdf", sep="", collapse=NULL),  width=6, height=6)
p + facet_wrap(~kraken2.db)
dev.off()

# ----- DIVERSITY -----

sample_data(fx0)$primer
sample_data(fx0)$ena.id

pAlpha = plot_richness(fx0,
                       shape = "library",
                       color = "primer",
                       measures = c("Observed", "Chao1", "Shannon", "InvSimpson"),
                       title = "OSD Alpha Diversity") +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, vjust = -0.001))
pdf(file = paste(img.path, today, "_fx0_diversityplot.pdf", sep="", collapse=NULL),  width=6, height=6)
pAlpha + geom_point(size = 5) 
dev.off()

# ------ SILVA ONLY -----

fx0.silva <- subset_samples(fx0, kraken2.db=="silva")

# collapse taxonomy to Kingdom rank
data.kingdoms <- tax_glom(fx0.silva, taxrank=rank_names(fx0.silva)[1], 
                          NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

p = plot_bar(data.kingdoms, x="library", fill = "Kingdom") +
  geom_bar( aes(color = Kingdom, fill = Kingdom ), 
            stat = 'identity', 
            position = 'stack') +
  theme_pubclean()

pdf(file = paste(img.path, today, "_fx0_silva.pdf", sep="", collapse=NULL),  
    width=6, height=6)
p
dev.off()

# ------ BACTERIA ----- 
bakt <- subset_taxa(fx0.silva, Kingdom=="Bacteria")
set_palette(p, "rickandmorty")
# collapse taxonomy to phylum rank
p.bakt <- tax_glom(bakt, 
                   taxrank=rank_names(bakt)[2], 
                   NArm=TRUE, 
                   bad_empty=c(NA, "", " ", "\t"))

p = plot_bar(p.bakt, x="library", fill = "Phylum") +
  geom_bar( aes(color = Phylum, fill = Phylum ), 
            stat = 'identity', 
            position = 'stack') +
  theme_pubclean()
pdf(file = paste(img.path, today, "_fx0_silva_bakt.pdf", sep="", collapse=NULL),  
    width=6, height=6)
p
dev.off()

# ------ EUKARYOTES ----- 

euk <- subset_taxa(fx0.silva, Kingdom=="Eukaryota")

# collapse taxonomy to phylum rank
p.euk <- tax_glom(euk, 
                   taxrank=rank_names(euk)[2], 
                   NArm=TRUE, 
                   bad_empty=c(NA, "", " ", "\t"))

p = plot_bar(p.euk, x="library", fill = "Phylum") +
  geom_bar( aes(color = Phylum, fill = Phylum ), 
            stat = 'identity', 
            position = 'stack') 

pdf(file = paste(img.path, today, "_fx0_silva_euk.pdf", sep="", collapse=NULL),  
    width=6, height=6)
p
dev.off()

# ----- DIVERSITY -----


pAlpha = plot_richness(fx0.silva,
                       shape = "library",
                       color = "stn.id",
                       measures = c("Observed", "Chao1", "Shannon", "InvSimpson"),
                       title = "OSD Alpha Diversity")
pdf(file = paste(img.path, today, "_fx0_silva_diversity.pdf", sep="", collapse=NULL),  
    width=6, height=6)
pAlpha + geom_point(size = 5)
dev.off()

# here it is meaningless to use beta-diversity
# because we compare the same sample!

# ----- RAREFACTION -----

fx0.silva.rare <- rarefy_even_depth(fx0.silva, replace=TRUE)

# collapse taxonomy to phylum rank
fx0.silva.rare.k <- tax_glom(fx0.silva.rare, 
                  taxrank=rank_names(fx0.silva.rare)[1], 
                  NArm=TRUE, 
                  bad_empty=c(NA, "", " ", "\t"))

p = plot_bar(fx0.silva.rare.k, x="primer", fill = "Kingdom") +
  geom_bar( aes(color = Kingdom, fill = Kingdom ), 
            stat = 'identity', 
            position = 'stack') 

pdf(file = paste(img.path, today, "_fx0_silva_rare_comp.pdf", sep="", collapse=NULL),  
    width=6, height=6)
p
dev.off()

# ----- BACTERIA RAREFACTION -----

bakt <- subset_taxa(fx0.silva, Kingdom == "Bacteria")
bakt1 <- subset_samples(bakt, ena.id=="ERR867776" | ena.id=="ERR771024")

bakt.rare <- rarefy_even_depth(bakt1, replace=TRUE)

# collapse taxonomy to phylum rank
p.bakt.rare <- tax_glom(bakt.rare, 
                   taxrank=rank_names(bakt.rare)[2], 
                   NArm=TRUE, 
                   bad_empty=c(NA, "", " ", "\t"))

p = plot_bar(p.bakt.rare, x="primer", fill = "Phylum") +
  geom_bar( aes(color = Phylum, fill = Phylum ), 
            stat = 'identity', 
            position = 'stack') 

pdf(file = paste(img.path, today, "_fx0_silva_rare_bakt.pdf", sep="", collapse=NULL),  
    width=6, height=6)
p
dev.off()

# ----- EUKARYOTES RAREFACTION -----

euk <- subset_taxa(fx0.silva, Kingdom == "Eukaryota")
euk1 <- subset_samples(euk, ena.id=="ERR867777" | ena.id=="ERR771024")

euk.rare <- rarefy_even_depth(euk1, replace=TRUE)

# collapse taxonomy to phylum rank
p.euk.rare <- tax_glom(euk.rare, 
                        taxrank=rank_names(euk.rare)[2], 
                        NArm=TRUE, 
                        bad_empty=c(NA, "", " ", "\t"))

p = plot_bar(p.euk.rare, x="primer", fill = "Phylum") +
  geom_bar( aes(color = Phylum, fill = Phylum ), 
            stat = 'identity', 
            position = 'stack') 


pdf(file = paste(img.path, today, "_fx0_silva_rare_euk.pdf", sep="", collapse=NULL),  
    width=6, height=6)
p
dev.off()

# ----- DIVERSITY -----


pAlpha = plot_richness(bakt.rare,
                       shape = "library",
                       color = "stn.id",
                       measures = c("Observed", "Chao1", "Shannon", "InvSimpson"),
                       title = "OSD Alpha Diversity")

pdf(file = paste(img.path, today, "_fx0_silva_rare_bakt_diversity.pdf", sep="", collapse=NULL),  
    width=6, height=6)
pAlpha + geom_point(size = 5)
dev.off()

pAlpha = plot_richness(euk.rare,
                       shape = "library",
                       color = "stn.id",
                       measures = c("Observed", "Chao1", "Shannon", "InvSimpson"),
                       title = "OSD Alpha Diversity")

pdf(file = paste(img.path, today, "_fx0_silva_rare_euk_diversity.pdf", sep="", collapse=NULL),  
    width=6, height=6)
pAlpha + geom_point(size = 5)
dev.off()
