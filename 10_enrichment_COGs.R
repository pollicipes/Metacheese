
library(data.table)
library(dplyr)
library(tidyr)

setwd("/Users/jrodriguez/Projects/metacheese/anvio_dirs/pangenome")

ax = c("lactis","cremoris")
tp = "FUNCTION" # FUNCTION
pv = 0.01
a = fread(paste0("enriched-COG20_",tp,"_subspecies.txt"))
a$id = seq(1,dim(a)[1])

# DEFINE A FUNCTION TO FIND THOSE POSITIONS WHERE THE STRAINS ARE PRESENT
# Returns a vector with indexes regarding the initial list

get_index <- function(gr = gr,a = a){
  expres = paste0("^",gr,"$")
  # h = grep(gr, a$associated_groups, value=TRUE, fixed = T)
  h = grep(gr, a$associated_groups, fixed = T)
  entries = a[c(h),]
  entries$id_sel = seq(1,dim(entries)[1])
  # dim(entries)
  assoc_groups = entries$associated_groups
  lh = strsplit(x = assoc_groups, split = ", ")
  result = lapply(lh, function(x) grep(expres, x))
  resx = lapply(1:length(result),
                 function(i)
                 {if(length(unlist(result[[i]])) != 0){
                   print(i)
                 }
                   })
  hindex = unlist(resx)
  sel_hindex = entries[hindex,]$id
  return(sel_hindex)
}

# USE THE FUNCTION TO GET THE INDEXES BACK
index.list = list()
for (gr in ax){
  x1 = get_index(gr=gr, a=a)
  index.list[[gr]] = x1
}

# Uniques in LACTIS
lind = setdiff(index.list$lactis,index.list$cremoris)
lactis_uniq = a %>%  
  slice(lind) %>%
  filter(adjusted_q_value < pv) %>%
  select(accession)

# Uniques in CREMORIS
cind = setdiff(index.list$cremoris,index.list$lactis)
cremoris_uniq = a %>%  
  slice(cind) %>%
  filter(adjusted_q_value < pv) %>%
  select(accession)

cogs = fread("/Users/jrodriguez/opt/miniconda3/envs/anvio-7.1/lib/python3.6/site-packages/anvio/data/misc/COG/COG20/COG.txt", header = FALSE)
colnames(cogs) = c("COG","group","description","pathway")

# SELECT THE GROUPS WHERE THEY BELONG ####

# Lactis
# View(lactis_uniq_x)
lactis_uniq_x = unlist(strsplit(lactis_uniq$accession,"!!!"))
lactis_cogs = cogs[cogs$COG %in% lactis_uniq_x,]
lactis_cogs_final = lactis_cogs %>%
  mutate(col1 = strsplit(as.character(group), ", ")) %>%
  unnest(col1) %>%
  filter(col1 != "") %>%
  select(COG, col1, description, pathway) %>%
  rename(group = col1)


# Cremoris
cremoris_uniq_x = unlist(strsplit(cremoris_uniq$accession,"!!!"))
cremoris_cogs = cogs[cogs$COG %in% cremoris_uniq_x,]
cremoris_cogs_final = cremoris_cogs %>%
  mutate(col1 = strsplit(as.character(group), ", ")) %>%
  unnest(col1) %>%
  filter(col1 != "") %>%
  select(COG, col1, description, pathway) %>%
  rename(group = col1)

# Lactis SELECTION 
# View(lactis_cogs_final)
tab_lactis = lactis_cogs_final %>% filter(group %in% c("P","I","K","H","G","C","E")) %>% print(n=300) 
write.table(tab_lactis, "~/Projects/metacheese/anvio_dirs/IMAGES/COGplot_pv_lactis.csv",sep="\t",quote = F)
lactis_cogs_final %>% filter(pathway != "") %>% select(pathway)
lactis_cogs_final %>% arrange(group,description) %>% print(n=200)

# Cremoris SELECTION
# View(cremoris_cogs_final)
# cremoris_cogs_final %>% filter(group == "L")
tab_cremoris = cremoris_cogs_final %>% filter(group %in% c("L","J","Q","V","F")) %>% print(n=300)
write.table(tab_cremoris, "~/Projects/metacheese/anvio_dirs/IMAGES/COGplot_pv_cremoris.csv",sep="\t",quote = F)
cremoris_cogs_final %>% arrange(group) %>% print(n=200)

categ = fread("/Users/jrodriguez/opt/miniconda3/envs/anvio-7.1/lib/python3.6/site-packages/anvio/data/misc/COG/COG20/CATEGORIES.txt", header = FALSE)
colnames(categ) = c("group","description")

lcounts_tmp = table(lactis_cogs_final$group)
lnames = names(lcounts_tmp)
missl = setdiff(toupper(letters),lnames)
zerol = rep(1,length(missl))
names(zerol) = missl
lcounts = sort(c(zerol,lcounts_tmp))
lcounts.sort = lcounts[order(names(lcounts))]

ccounts_tmp = table(cremoris_cogs_final$group)
cnames = names(ccounts_tmp)
missc = setdiff(toupper(letters),cnames)
zeroc = rep(1,length(missc))
names(zeroc) = missc
ccounts = sort(c(zeroc,ccounts_tmp))
ccounts.sort = ccounts[order(names(ccounts))]

# FOLD CHANGE
log2FC = log2(lcounts.sort) - log2(ccounts.sort)
# Decreasing sort
fc = sort(log2FC,decreasing = TRUE)
fc_z = reshape2::melt(fc)
fc_z$names = rownames(fc_z)
fc_z$names = factor(fc_z$names,levels = rev(rownames(fc_z)))
categ = categ[match(rownames(fc_z), rev(categ$group)),]

fc_x = merge(fc_z,categ, by.x = "names", by.y = "group")
fc_x = fc_x[order(fc_x$value),]

# PLOT
library(ggplot2)
fun_color_range = colorRampPalette(c("#1b98e0", "firebrick2"))
cols = fun_color_range(26)

fc_x$description = factor(fc_x$description,
                levels = fc_x$description)

cogplot = ggplot(data = fc_x, aes(x = description, y = value, fill = names)) +
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "grey", size=0.5) +
  geom_bar(stat = "identity") + 
  theme_minimal() + 
  xlab("COG Categories") + 
  ylab("Significant log2FC lactis:cremoris") + 
  scale_fill_manual(values = cols) +
  # scale_x_discrete(labels=fc_x$description)  +
  coord_flip() + 
  theme(legend.position = "none") +
  ggtitle("Enrichment of COG functions L. lactis vs L. cremoris")

cogplot
ggsave(filename = paste0("~/Projects/metacheese/anvio_dirs/IMAGES/COGplot_pv",pv,".pdf"),plot = cogplot,device = "pdf")

# FC = 2 ^ log2FC







