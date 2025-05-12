#### FOR THE ANOVA ANALYSES: ABUNDANCE LACTOCOCCUS - RIND + CORE

# !!! USE R 4.x.x version # 

library(tidyverse);
library(ggpubr);
library(dplyr);
library(ggplot2);
library(car);
library(pwr2);
library(data.table)
library(sjstats)
library(effectsize)

# LOAD DATA ####
setwd(paste0("~/HPC_data/anvio_dirs/runs/"))
norm_type = "proportion"
load(file = paste0("~/Projects/metacheese/anvio_dirs/df.tx_object_abundances_norm_",norm_type,".Rda")) ## Object df.tx, with the % abundances per each MAG
# df.tx
# index order in data frame
ord = c(14,13,12,11,2,1,6,5,4,3,10,9,8,7)
# labels for the data frame
xlbs = c("preinoculation", "postinoculation", "beforeMolding", "beforeBrining", "T0_10C_Rind", "T0_10C_Core", "T1_7d_13C_Rind", "T1_7d_13C_Core", "T1_35d_13C_Rind", "T1_35d_13C_Core", "T2_7d_20C_Rind", "T2_7d_20C_Core", "T2_35d_20C_Rind", "T2_35d_20C_Core")

# Select the corresponding entries
#### samps = c("Clostridium_tyrobutyricum")
samps = c("LAC2","CRE6")
all_bts = data.frame()
for (bt in c("D11","D12","D13")){
  selx = df.tx %>% as.data.frame %>% select(contains(bt)) %>% select(all_of(ord)) %>% 
    filter(rownames(df.tx) %in% samps) %>% 
    transpose %>% 
    mutate(batch = rep(bt, 14), step=xlbs)
  all_bts = bind_rows(all_bts,selx)
}
core = all_bts %>% filter(grepl("Core", step)) %>% mutate(lacto = V1+V2) %>% select(3,4,5) %>% mutate(site = "Core")
rind = all_bts %>% filter(grepl("Rind", step)) %>% mutate(lacto = V1+V2) %>% select(3,4,5) %>% mutate(site = "Rind")

my_data = bind_rows(core,rind)
my_data$site = factor(my_data$site, levels = c("Core","Rind"))
my_data$batch = factor(my_data$batch, levels = c("D11","D12","D13"))

# ANOVA TWO-WAY ####

dx = ggboxplot(my_data, x = "batch", y = "lacto", color = "site",
          palette = c("#00AFBB", "#E7B800")) + 
  xlab("Batch") + 
  ylab("Total Lactococcus abundance %")
dx
ggsave("~/Projects/metacheese/anvio_dirs/IMAGES/SUPPL_FigX_Rind-Core_abundance.pdf",dx)

## IMPACT OF EACH FACTOR ####
res.aov2 = aov(lacto ~ site + batch, data = my_data)
summary(res.aov2)
effectsize::omega_squared(res.aov2, partial = FALSE)

### POWER ANALYSIS ####
# 3 groups for batch; 2 groups for site. 10 elements per batch, 15 elements per site.
pwr.2way(a=3, b=2, alpha=0.05, size.A=10, size.B=15, 
         f.A=0.71, f.B=0.07)

## INTERACTION BETWEEN FACTORS ####
res.aov3 = aov(lacto ~ site * batch, data = my_data)
summary(res.aov3)
effectsize::omega_squared(res.aov3, partial = FALSE)

# ASSUMPTIONS OF THE MODEL ####
# Residual distributions
plot(res.aov3, 1)
# Homogeneity of variances.
leveneTest(lacto ~ site*batch, data = my_data)
# Normality assumption
plot(res.aov3, 2)
aov_residuals <- residuals(object = res.aov3)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )

# SUMMARY TABLE WITH THE MEAN AND SE PER GROUP 
group_by(my_data, site, batch) %>%
  summarise(
    count = n(),
    mean = mean(lacto, na.rm = TRUE),
    sd = sd(lacto, na.rm = TRUE)
  )

# Interaction plot
interaction.plot(x.factor = my_data$batch, trace.factor = my_data$site, 
                 response = my_data$lacto, fun = mean, 
                 type = "b", legend = TRUE, 
                 xlab = "Batch", ylab="Lactococcus abundance",
                 pch=c(1,19), col = c("#00AFBB", "#E7B800"))

# ANOVA 2 WAY TEMPERATURE EFFECT.

# ANOVA 2 WAY Clostridium.






