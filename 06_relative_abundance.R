
library(tidyverse)
library(readxl)
library(phyloseq)
library(ape)
library(data.table)
library(microbiome)
library(RColorBrewer)
library(ggpubr)
library(scales)
# library(rstatix)
options(scipen = 10)

# samp="jacob_example"
# setwd(paste0("~/HPC_data/anvio_dirs/runs/",samp))
# # FOR OFFLINE TESTING:
# setwd(paste0("~/Downloads/for_MAG_testing/"))

### MY DATA ####
setwd(paste0("~/HPC_data/anvio_dirs/runs/"))

# COVERAGE
df0 = read.table("MAG_coverage_ALL.txt", header=T)
# df0 = read.table("D11/07_RELATIVE_ABUNDANCE/SUMMARIZED/bins_across_samples/mean_coverage.txt", header=T)
read_cvg_table <- function(df0 = df0) {
  df = as.matrix(df0)
  rownames(df) = df[,1]
  df = df[,-c(1)]
  class(df) = "numeric"
  return(df)
  
}
df = read_cvg_table(df0)

# BINS
bin_sum = fread("bins_summary.txt",header=T)
# genome size
gsize = data.frame(MAG = bin_sum$bins, Length = bin_sum$total_length)

# TAXONOMY
tax = fread("MAG_Tax_ALL.txt")
tax = as.data.frame(tax)
# Tidy up tax
tax = as.matrix(tax)
rownames(tax) = rownames(df)
tax = tax[,-c(1)]

# METADATA
md = fread("metadata.txt")
# Tidy up md
md = as.matrix(md)
rownames(md) = md[,1]
md = as.data.frame(md)
# View(md)

# GET NUMBER OF READS PER MAG/SAMPLE
read_count <- function(df = df, gsize = gsize, rleng = NULL){
  if(is.null(rleng)){
    stop("Define the read length.")
  }
  ln = as.matrix(gsize$Length)
  dfx = sweep(df, 1, STATS = ln, FUN = "*")
  dfx1 = round(dfx/rleng) # Divide by read length
  return(dfx1)
  }
# This is a Lander-Waterman transformation to get read count out from coverage. https://www.biostars.org/p/107070/
df.x = read_count(df,gsize,rleng = 100)

# Export numbers to excel format
# write_excel_csv2(as.data.frame(df.x), "~/Projects/metacheese/drafts/FINAL_FIGURES/Table_read_counts.csv")

# READS RECRUITED BY OUR MAGS
# Adding the proportion of reads mapped, out of the total metagenome.
tot_metag = read.table("all_read_counts.txt")
colnames(tot_metag) = c("step","batch","reads","path")
tot_metag$id = paste0(tot_metag$step,"_",tot_metag$batch)
sumas = as.data.frame(colSums(df.x))
sumas$id = rownames(sumas)
colnames(sumas)=c("mapped","id")
prop_recruit1 = merge(tot_metag,sumas, by = c("id"))
prop_recruit1$perc = (prop_recruit1$mapped / prop_recruit1$reads)*100

# RATIO OF VIRUS TO HOST (BACTERIA)
vir_prop <- function(x){
  vir_cov = sum(x[names(x) %in% c("Ceduovirus_1","Ceduovirus_2","Skunavirus_2")])
  lacto_cov = sum(x[names(x) %in% c("LAC2","CRE6")])
  ratio = vir_cov/lacto_cov
  return(ratio)
}
vir_tab = as.data.frame(apply(df, 2, vir_prop))
vir_tab$id = rownames(vir_tab)
colnames(vir_tab) = c("phages_per_cell","id")
prop_recruit2 = merge(prop_recruit1, vir_tab, by.x = "id" )
# View(prop_recruit2)

# STATS ON VIRAL PROPORTIONS AND READS ####
st_exc = c("milk", "preinoculation","postinoculation","beforeBrining","beforeMolding")
# Percentage of reads mapped back to MAGs:
mean(prop_recruit2[!(prop_recruit2$step %in% st_exc),]$perc)

# Phages per cell: DK vs. SWE
dk = c("D11","D12","D13")
dk.p = prop_recruit2[prop_recruit2$batch %in% dk,]
dk.ph = dk.p[!(dk.p$step %in% st_exc),]$phages_per_cell
## Only if we were to include the other Sweden plants.
swe = c("Kalmar93","Kalmar94")
# swe.p = prop_recruit2[prop_recruit2$batch %in% swe,]
# swe.ph = swe.p[!(swe.p$step %in% st_exc),]$phages_per_cell
# wilcox.test(swe.ph,dk.ph)
# mean(dk.ph)/mean(swe.ph)

# PHAges ON D11 vs D12, D13
prop_recruit2
# Proportion recruited by the phages in the microfiltrated batch.
phag_11 = prop_recruit2[prop_recruit2$batch %in% "D11",]
# Proportion recruited by the phages in the other ones. 
phag_x = prop_recruit2[prop_recruit2$batch %in% c("D12","D13"),]
p11 = phag_11[!(phag_11$step %in% st_exc),]$phages_per_cell
px =  phag_x[!(phag_x$step %in% st_exc),]$phages_per_cell
mean(px)/mean(p11)
message("Proportion of phages reads recruited on average by the non-microfiltrated is > 2.2x higher than in the microfiltration trial.\n")

# SELECTING OTHER NON-MAIN MAGS: ####
# Removes phages and Lactococcus MAGs
phages_rm = TRUE
if(isTRUE(phages_rm)){
  df.x = df.x[-c(1,2,3,6),]
}

# NORMALIZATION ####
# TPM normalisation
NormalizeTPM <- function(data, sce = NULL, tr_length = NULL, log = FALSE, scale = 1e+06, pseudo.count = 1) {
  
  # Median length of all transcripts for a given gene
  med_length <- stats::aggregate(x = tr_length$transcript_length,
                                 by = list(hgnc_symbol = tr_length$hgnc_symbol), FUN = stats::median)
  common <- intersect(rownames(data), med_length$hgnc_symbol)
  if (length(common) < 2)
    stop("No enough overlap between data rownames and transcript length.\n")
  data <- data[common, ]
  
  med_length <- med_length[match(common, med_length$hgnc_symbol), ]
  
  # divide by length of transcript in kb - reads per kilobase (RPK)
  data1 <- sweep(data, 1, STATS = med_length$x/1000, FUN = "/")
  # divide by million RPK in sample - transcripts per million (TPM)
  data2 <- sweep(data1, 2, STATS = colSums(data1)/(10^6), FUN = "/")
  data2[,1]
  data3 <- data2/scale
  
  if (log) {
    if (pseudo.count == 0)
      warning("Using 0 pseudocount: Inf may be generated.\n")
    data <- log2(data + pseudo.count)
  }
  
  if(is.null(sce)){ return(data3) } else{ return(sce) }
}

## DEFINE BATCH AND NORM TYPE ##
n_batch = "D13"
# norm_type = "rawcounts"
# norm_type = "proportion"
# norm_type = "TPM"
for (norm_type in c("rawcounts","proportion","TPM")){
  if (norm_type == "proportion") {
    # Alternative 1 ####
    # We can transform the reads counts into compositional abundance (all sum equal) (Lander-Waterman)
    df.tx = microbiome::transform(df.x, transform = "compositional")
    # save(df.tx, file = paste0("~/Projects/metacheese/anvio_dirs/df.tx_object_abundances_norm_",norm_type,".Rda"))
    tag = "Relative read proportion"
    } else if (norm_type == "rawcounts") {
        # using raw read counts:
        df.tx = df.x
        tag = "Raw read count"
    } else if (norm_type == "TPM"){
      # Alternative 2 ####
        size = data.frame(MAG = bin_sum$bins, Length = bin_sum$total_length)
        colnames(size) = c("hgnc_symbol", "transcript_length")
        df.tx = NormalizeTPM(df.x, tr_length = size, scale = 1e+06)
        tag = "Normalized relative abundance"
    }
  
  # PHYSEQ OBJECTS ####
  physeq.tmp = phyloseq(otu_table(df.tx, taxa_are_rows=TRUE),
                        tax_table(tax),
                        sample_data(md))
  random_tree = rtree(ntaxa(physeq.tmp), rooted=TRUE, tip.label=taxa_names(physeq.tmp))
  physeq.tmp = merge_phyloseq(physeq.tmp, random_tree)
  
  # PLOTS ####
  ##### Subset batch #####
  steps = read.table("steps.txt",header=T)
  n_batch_levs = as.character(steps[,c(n_batch)],na.rm=T)
  n_batch_levs = n_batch_levs[!is.na(n_batch_levs)]
  MAGs_ps = subset_samples(physeq.tmp, batch==n_batch)
  # reordering levels
  sample_data(MAGs_ps)$ID = factor(sample_data(MAGs_ps)$ID, levels = n_batch_levs)
  # renaming levels
  xlbs = gsub(paste0("_",n_batch),"",n_batch_levs) # removing the name out of the step
  levels(sample_data(MAGs_ps)$ID) = xlbs
  sample_data(MAGs_ps)$place = factor(sample_data(MAGs_ps)$place, levels = c("General","Rind","Core"))
  sample_data(MAGs_ps)$temperature = factor(sample_data(MAGs_ps)$temperature, levels = c("RT","10C","13C","20C"))
  set.seed(1) # This makes the plot appear the same each time it is run 
  # Colors
  my_pal_init = RColorBrewer::brewer.pal(n = 11, name = "Set3")
  # Select to be consistent with previous color codes
  my_pal = my_pal_init[-c(6,10)]
  # Plot adjustments if we don't want the phages/lactococcus
  if(isTRUE(phages_rm)){
    # scales::show_col(my_pal)
    my_pal = my_pal_init[-c(1,2,7,8)]
  }
  
  #### PLOTTING ####
  p_gut = plot_bar(MAGs_ps, x = "ID", fill = "t_species") + 
    geom_bar(aes(color=t_species, fill=t_species), stat="identity", position="stack") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, size=15,hjust=1,vjust=0.3),plot.margin = margin(0, 0, 0, 0, "cm"),
          axis.text.y = element_text(size=15),
          axis.title.y = element_text(size=20),
          axis.title.x = element_text(size=20)) +
    
    # scale_x_discrete(labels = xlbs) +
    scale_fill_manual(values = my_pal) +
    scale_color_manual(values = my_pal) +
    labs(x = "Time point") +
    # xlab(NULL) +
    labs(y = tag )
  p_gut 
  
  # SAVING PLOTS ####
  if(isTRUE(phages_rm)){
    fname = paste0("~/HPC_data/anvio_dirs/runs/IMAGES_Abundance/PHAGES+LACTOC_REMOVED_",n_batch,"_",tag,".pdf")
  } else {
    fname = paste0("~/HPC_data/anvio_dirs/runs/IMAGES_Abundance/",n_batch,"_",tag,".pdf")
  }
  if(norm_type == "rawcounts") {
    p = p_gut + facet_wrap(~temperature, scales = "free", nrow = 1) + 
      theme(legend.position="bottom",strip.text.x = element_blank()) + 
      scale_y_continuous(limits = c(0,max(colSums(MAGs_ps@otu_table)) + 100000), labels = scientific) + 
      labs(color="Species",fill="Species") #,title = n_batch, subtitle = tag)
    if(isTRUE(phages_rm)){
      ggsave(filename = fname, plot = p)
      next;
      }
  } else {
    p = p_gut + facet_wrap(~temperature, scales = "free", nrow = 1) + 
      theme(legend.position="bottom",
            legend.text = element_text(size = 12),strip.text.x = element_blank()) + 
      labs(color="Species",fill="Species") #,title = n_batch, subtitle = tag)  
    if(isTRUE(phages_rm)){
      ggsave(filename = fname, plot = p)
      next;
      }
    # ggsave(filename = fname, plot = p)
  }
  
  ##### PLOTTING READS RECRUITED & VIRUS #####
  ###### Reads ####
  prop_recruit = merge(prop_recruit2[prop_recruit2$batch == n_batch,], md, by.x = "id",by.y = "ID")
  prop_recruit$id = factor(prop_recruit$id, levels=n_batch_levs)
  prop_recruit$temperature = factor(prop_recruit$temperature, levels=c("RT","10C","13C","20C"))
  
  p_reads = ggplot(data = prop_recruit, aes(x=id, y=perc)) +
    geom_bar(stat="identity", fill="gray")+
    theme_classic() +
    theme(axis.title.y = element_text(size=20), 
          strip.text.x = element_text(size = 17),
          axis.text.y = element_text(size=15)) +
    # theme(axis.text.x = element_blank()) +
    scale_y_continuous(expand = expansion(mult = c(0, .005)), limits = c(0, 100)) +
    xlab("") +
    ylab("% reads recruited\n over metagenome") + 
    facet_wrap(~temperature, scales = "free", nrow = 1) +
    labs(title = n_batch) + # , # subtitle = "Reads recruited to MAGs") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  p_reads
  fname_r = paste0("~/HPC_data/anvio_dirs/runs/IMAGES_Abundance/",n_batch,"_reads_recruited.pdf")
  # ggsave(filename = fname_r, plot = p_reads)
  
  ###### Virus ####
  p_vir = ggplot(data = prop_recruit, aes(x=id, y=log2(phages_per_cell))) +
    geom_bar(stat="identity", fill="lightblue")+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust=1), 
          axis.title.y = element_text(size=20),
          axis.text.y = element_text(size=15)) +
    # theme(axis.text.x = element_blank()) +
    # ylim(0,100) + 
    xlab("") +
    ylab("Log2 ratio \n phages : Lactococcus") + 
    facet_wrap(~temperature, scales = "free", nrow = 1) +
    theme(strip.background = element_blank(), strip.text.x = element_blank()) +
    # scale_y_continuous(limits = c(0,max(prop_recruit2$phages_per_cell))) +
    scale_y_continuous(limits = c(-8,8), expand = c(0,0)) +
    # labs(title = n_batch, subtitle = "Phage copies per Lactococcus cell") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    geom_hline(yintercept=0,linetype="dashed",color="darkred",size=0.35)
  
  # p_vir
  fname_v = paste0("~/HPC_data/anvio_dirs/runs/IMAGES_Abundance/",n_batch,"_phages_copies.pdf")
  # ggsave(filename = fname_v, plot = p_vir)
  
  # Save plot
  full = ggarrange(p_reads, p_vir, p, ncol = 1, nrow = 3, align = "v", heights = c(0.12,0.12,0.7))
  # ggsave(filename = paste0("~/HPC_data/anvio_dirs/runs/IMAGES_Abundance/FULL ",n_batch," ",tag,".pdf"), plot = full)
  ggsave(filename = paste0("~/Projects/metacheese/drafts/IMAGES_SUPPL_TABLES_Metacheese_190723/FULL ",n_batch," ",tag,"_EDIT_VERSION.pdf"), plot = full)
  }

stop("DONE")

# RELATION BETWEEN CREMORIS AND LACTIS #### 

library(dplyr)

# Column order
ord = c(14,13,12,11,2,1,6,5,4,3,10,9,8,7)
# Step Names
print(xlbs)
# Samples to use
samps = c("LAC2","CRE6") 

# Function to extract abundances 
extract_abundance <- function(df.tx = df.tx, ord = ord, samps = samps, xlbs = xlbs){
  # Emtpty df to store abundances
  all_bts=data.frame()
  for (bt in c("D11","D12","D13")){
    selx = df.tx %>% as.data.frame %>% select(contains(bt)) %>% select(ord) %>% 
      filter(rownames(df.tx) %in% samps) %>% 
        transpose %>% 
        mutate(batch = rep(bt, 14), step=xlbs)
    all_bts = bind_rows(all_bts,selx)
  }
  return(all_bts)
}
all_bts = extract_abundance(df.tx = df.tx, ord = ord, samps = samps, xlbs = xlbs)

# Proportion lactis/cremoris
colnames(all_bts) = c("cremoris","lactis","batch","step")
all_bts$log2FC = log2(all_bts$lactis) - log2(all_bts$cremoris)

# Reorder levels
all_bts$step = factor(all_bts$step,levels=xlbs)

# PLOTTING ####

# Lactis / cremoris #####

pit = ggplot(data = all_bts,aes(x=step,y=log2FC,colour = batch)) +
  geom_hline(yintercept = 0, color="gray",lty="dashed") +
  geom_point(shape = 19,size=2.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size=10,hjust=1,vjust=0.3)) +
  scale_colour_brewer(palette = "Set2") +
  labs(x = "Time point", y = "log2FC lactis : cremoris") +
  ggtitle("Proportion of lactis vs. cremoris abundances")
pit
ggsave("~/Projects/metacheese/anvio_dirs/IMAGES/FigX_proportion_lactis_cremoris.pdf")  
  

# Virus/strain Lactococcus #####

samps = c("Ceduovirus_1","Ceduovirus_2","LAC2","CRE6")
all_bts = extract_abundance(df.tx = df.tx, ord = ord, samps = samps, xlbs = xlbs)
colnames(all_bts) = c("cremoris","Ceduovirus_1","Ceduovirus_2","lactis","batch","step")

# Calculate proportions
all_bts$fc_lact = log((all_bts$lactis + 0.01)/(all_bts$Ceduovirus_1 + 0.01))
all_bts$fc_crem = log((all_bts$cremoris + 0.01)/(all_bts$Ceduovirus_2 + 0.01))
# all_bts$fc_lact = log2(all_bts$lactis + 0.01) - log2(all_bts$Ceduovirus_1 + 0.01)
# all_bts$fc_crem = log2(all_bts$cremoris + 0.01) - log2(all_bts$Ceduovirus_2 + 0.01)

# Reorder levels
all_bts$step = factor(all_bts$step,levels=xlbs)
all_bts$period = rep(c(rep("Starting",4),rep("Temp0",2),rep("Ripening",8)),3)

# Correlation #####
truc2 = ggplot(all_bts, aes(x = fc_lact, y = fc_crem, shape = batch, colour = period)) +
  geom_point(size=c(2)) +
  scale_shape_manual(values=c(3, 13, 17)) +
  theme_classic() +
  scale_colour_brewer(palette = "Set1") +
  xlim(-0.2,5) + 
  ylim(-0.2,5) +
  coord_equal() +
  labs(x = "L.lactis / Ceduovirus_1", y = "L.cremoris / Ceduovirus_2") +
  ggtitle("Correlation Lactococcus strains/virus") +
  geom_smooth(method=lm,aes(group=period),se=FALSE) +
  ggpubr::stat_cor(method = "pearson",aes(group=period))
ggsave("~/Projects/metacheese/anvio_dirs/IMAGES/SUPPL_FigX_Correlation_bacteriaStrain_phage.pdf",plot = truc2)  

# Proportion #####

# Adapt dataframe to plot
# all_fc = reshape2::melt(all_bts[,c("batch","step","fc_lact","fc_crem")])
# colnames(all_fc) = c("batch","step","strain","fold")

# truc = ggplot(data = all_fc, aes(x = step, y = fold, colour = strain, shape = batch)) +
#   geom_hline(yintercept = 1, color="gray",lty="dashed") +
#   geom_point() + # shape=c(19,17),size=2.5) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, size=10,hjust=1,vjust=0.3)) +
#   scale_colour_brewer(palette = "Set2") +
#   labs(x = "Time point", y = "Ratio Lactococcus/virus") +
#   ggtitle("Ratio  Lactococcus strain to phage")
# truc

# Average phages per batch (through ripening) ####
prop_recruit2 %>% filter(!step %in% c("preinoculation","postinoculation","beforeMolding","beforeBrining")) %>% group_by(batch) %>% dplyr::summarize(Mean = mean(log2, na.rm=TRUE))

# Leuconostoc abundance ####

df.tx %>% as.data.frame %>% filter(rownames(df.tx) %in% "Leuconostoc_pseudomesenteroides") %>% select(contains("D13")) %>% transpose
# filter(!step %in% c("preinoculation","postinoculation","beforeMolding","beforeBrining")) %>% group_by(batch) %>% dplyr::summarize(Mean = mean(log2, na.rm=TRUE))

# END ####

stop("DONE UNTIL HERE")

####
#### It's correct until here 
####

#### OTHER STATISTICAL ANALYSES:
GP.ord <- ordinate(MAGs_ps, "PCoA", "unifrac", weighted=TRUE)
pcoa = plot_ordination(MAGs_ps, GP.ord, type="samples",
                       color="Location",
                       label="ID",
                       title="Microbiota composition of wild Norwegian salmon") +
  geom_point(size=5) + theme_bw() +
  scale_fill_manual(values = my_pal) +
  scale_color_manual(values = my_pal)
pcoa + facet_wrap(~Location, nrow = 1)

## Make physeq object
class(df) <- "integer"
df.rich <- round(df.tmp*1000000,0) # make TMP normalised data to round percentage
physeq.rich <- phyloseq(otu_table(df.rich,taxa_are_rows=TRUE),
                        tax_table(tax),
                        sample_data(md))
physeq.tmp
random_tree = rtree(ntaxa(physeq.tmp), rooted=TRUE, tip.label=taxa_names(physeq.tmp))
physeq.rich = merge_phyloseq(physeq.tmp, random_tree)
MAGs_div = subset_taxa(physeq.rich)


plot_richness(MAGs_div, x="Location", color="Location", measures = c("Observed","Shannon", "Simpson")) + 
  geom_point(size=5) + 
  xlab("Location") +
  scale_fill_manual(values = my_pal) +
  scale_color_manual(values = my_pal) + xlab("") +
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))
?plot_richness
richness <- estimate_richness(MAGs_div, measures = c("Observed","Shannon", "Simpson"))
md.richness <- md[match(rownames(richness),md$ID),]
richness <- cbind(md.richness, richness)
###



# TEST AREA ####


# gsize$Length * (/a)/

suma = sum(df.x[,1])
covg = df[,1]

a = (covg * gsize$Length) / (100 * suma)

data("esophagus")
View(esophagus@otu_table)
estimate_richness(esophagus, measures=c("Observed", "InvSimpson", "Shannon", "Chao1"))


data(dietswap)
x <- dietswap
x@otu_table[c(1:5),c(1:10)]
xt <- microbiome::transform(x, 'compositional')  
xt@otu_table[c(1:5),c(1:10)]
tot = sum(x@otu_table[,1])
head(x@otu_table[,1]/tot)

rr = sweep(df, 1, ln, `*`)
reads_mags = as.data.frame(rr/100)
reads_mags$sums = rowSums(reads_mags)
reads_mags$preinoculation/reads_mags$sums
mg1.t = sum((ln * mg1)/100)

# END ####