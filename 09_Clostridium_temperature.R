#### Abundance analyses; temperature effect on CLostridium.

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
# setwd(paste0("~/HPC_data/metacheese/anvio_dirs/runs/"))
norm_type = "proportion"
load(file = paste0("~/Projects/metacheese/anvio_dirs/df.tx_object_abundances_norm_",norm_type,".Rda")) ## Object df.tx, with the % abundances per each MAG
# index order in data frame
ord = c(14,13,12,11,2,1,6,5,4,3,10,9,8,7)
# labels for the data frame
xlbs = c("preinoculation", "postinoculation", "beforeMolding", "beforeBrining", "T0_10C_Rind", "T0_10C_Core", "T1_7d_13C_Rind", "T1_7d_13C_Core", "T1_35d_13C_Rind", "T1_35d_13C_Core", "T2_7d_20C_Rind", "T2_7d_20C_Core", "T2_35d_20C_Rind", "T2_35d_20C_Core")

# Select the corresponding entries
samps = c("Clostridium_tyrobutyricum")
# samps = c("LAC2","CRE6")
all_bts = data.frame()
for (bt in c("D11","D12","D13")){
  selx = df.tx %>% as.data.frame %>% select(contains(bt)) %>% select(all_of(ord)) %>% 
    filter(rownames(df.tx) %in% samps) %>% 
    transpose %>% 
    mutate(batch = rep(bt, 14), step=xlbs)
  all_bts = bind_rows(all_bts,selx)
}

all_bts = subset(all_bts, select = -c(batch))
all_bts$trial = c(rep("PT1",14), rep("PT2",14), rep("PT3",14))

dtx0 = all_bts %>% filter(trial %in% c("PT1","PT2","PT3")) %>% slice(7:14,21:28,35:42)
dtx0$site = c("Rind","Core")
dtx0$temp = c("13C","13C","13C","13C","20C","20C","20C","20C")
dtx0$days = c("7d","7d","35d","35d")

dtx0$site = factor(dtx0$site, levels = c("Rind","Core"))
dtx0$trial = factor(dtx0$trial, levels = c("PT1","PT2","PT3"))
dtx0$temp = factor(dtx0$temp, levels = c("13C","20C"))
dtx0$days = factor(dtx0$days, levels = c("7d","35d"))
dtx0$perc = dtx0$V1*100

dtx0 %>%
  group_by(trial) %>%
  summarize(clos_perc = sum(perc))


# ANOVA TWO-WAY ####
comparisons <- list(c("PT1", "PT2"), c("PT2", "PT3"), c("PT1", "PT3"))
tit = expression(paste(italic("Clostridium tyrobutyricum"), " abundance %"))

dx = ggboxplot(dtx0, x = "trial", y = "perc", fill = "trial") + 
  theme(text = element_text(size = 15),axis.title = element_text(face = 'italic')) +
  
  scale_fill_manual(values = c("tomato1","darkseagreen","dodgerblue")) +
  labs(y = tit) +
  xlab("Trial") + 
  # ggtitle(tit) +
  stat_compare_means(method="wilcox.test",comparisons = comparisons) +
  # stat_compare_means(label="p.format", comparisons = comparisons, method = "wilcox.test") +
  guides(color = "none")
  # stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),comparisons = comparisons)
dx
dev.off()
ggsave("~/Projects/metacheese/drafts/FINAL_FIGURES/SUPPL_FigX_Clostridium_abundance.pdf",dx)
ggsave("~/Documents/Talks&Conferences/MIFFI_Conference_09042025/SUPPL_FigX_Clostridium_abundance.pdf",dx)

mean(dtx0[dtx0$trial == "PT3",]$perc)/mean(dtx0[dtx0$trial == "PT1",]$perc)
mean(dtx0[dtx0$trial == "PT2",]$perc)/mean(dtx0[dtx0$trial == "PT1",]$perc)

# REMOVE THE PT1 AS THERE IS NO CLOSTRIDIUM
dtf = dtx0[dtx0$trial != "PT1",]

x7 = dtf[dtf$days == "7d",]$perc
x35 = dtf[dtf$days == "35d",]$perc
wilcox.test(x7, x35)

# PLOT BY DAYS
ggboxplot(dtf, x = "days", y = "perc", color = "days") + 
  theme(text = element_text(size = 16)) +
  # scale_color_manual(values = c("tomato1","darkseagreen","dodgerblue")) +
  xlab("Trial") + 
  stat_compare_means(method="wilcox.test") +
  ylab("Total C. tyrobutyricum abundance %")

# ANOVA FOR EACH FACTOR ####

# Temperature
res.aov2 = aov(perc ~ temp, data = dtf)
summary(res.aov2)
# effectsize::omega_squared(res.aov2, partial = FALSE)

# Days
res.aov2 = aov(perc ~ days, data = dtf)
summary(res.aov2)

# Clostridium abundance per temperature and days ####
days_labels <- c(`35d` = "Ripening 35 Days",`7d` = "Ripening 7 Days")
pt_temp = ggplot(dtf, aes(temp, perc)) +
  geom_boxplot(aes(fill = temp)) +
  facet_wrap(~days,labeller = as_labeller(days_labels)) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  stat_compare_means(method="wilcox.test") +
  guides(color="none") +
  theme_classic() +
  ylim(0,0.2) +
  theme(text = element_text(size = 15)) +
  # ggtitle(tit,subtitle = "Ripening time & temperature") +
  scale_x_discrete(labels = c("13 ºC","20 ºC")) + 
  ylab(tit) +
  xlab("Ripening temperature")
pt_temp
ggsave("~/Projects/metacheese/drafts/FINAL_FIGURES/SUPPL_FigX_Clostridium_RIPENING_DAYS.pdf",pt_temp)
ggsave("~/Documents/Talks&Conferences/MIFFI_Conference_09042025/Temperature_Clostridium_abundance.pdf",pt_temp)
# END ####