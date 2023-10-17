
###################################################################################

### Present results from comparisons in cox models 

###################################################################################

### FDR adjustment of associations for proteins in 20k

screen

R

library(ggplot2)
library(tidyverse)

# read in results 20k proteins 
comb <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/01_Cox_20k_results/basic.csv")
comb_full <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/01_Cox_20k_results/full.csv")
names(comb_full)[11] <- "max_tte"

# Do FDR correction
comb <- comb[order(comb$`P.Value`),]
comb$FDR <- p.adjust(comb$P.Value, method = "BH")

write.csv(comb, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/01_Cox_20k_results/basic_FDR.csv", row.names = F)

# Keep any that meet significance threshold 
keep1 <- comb[which(comb$FDR < 0.05),] 
keep = paste(keep1$Predictor, keep1$Outcome, sep = "_")

# Keep only those passing the threshold from the basic model in the fully adjusted 
comb_full$retain <- paste(comb_full$Predictor, comb_full$Outcome, sep = "_")
comb_full$col <- ifelse(comb_full$retain %in% keep & comb_full$P.Value < 0.05, "red", "black")

# Save a supplementary table with basic and full model results 
comb$retain <- paste(comb$Predictor, comb$Outcome, sep = "_")
join <- left_join(comb, comb_full, by = "retain")
join <- join[order(join$P.Value.x),]
join <- join[c("Predictor.x", "Outcome.x", "No..of.Cases.x", "No..of.Controls.x", "tte_mean_sd.x", 'max_tte.x',
  "Hazard.Ratio.x", "LCI.x", "UCI.x", "P.Value.x", "FDR", "cox.zph.x",
  "No..of.Cases.y", "No..of.Controls.y", "tte_mean_sd.y", "max_tte.y",
  "Hazard.Ratio.y", "LCI.y", "UCI.y",
  "P.Value.y", "cox.zph.y", "col")]
names(join) <- c("Marker", "Outcome", "Basic N Cases", "Basic N Controls", "Basic mean tte (sd)", "Basic max tte",
  "Basic HR", "Basic LCI", "Basic UCI", "Basic P", "Basic FDR P", "Basic zph",
  "Full N Cases", "Full N Controls", "Full mean tte (sd)", "Full max tte",
  "Full HR", "Full LCI", "Full UCI", "Full P", "Full zph",
  "Association")
write.csv(join, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/01_Cox_20k_results/Suppl_table_results.csv", row.names = F)


###################################################################################

### Plot HR comparison for both proteins

x <- comb_full

library(stringr)

x$Outcome <- str_replace(x$Outcome, "IHD", "Ischaemic Heart Disease")
x$Outcome <- str_replace(x$Outcome, "Diabetes", "Type 2 Diabetes")
x$Outcome <- str_replace(x$Outcome, "Stroke", "Ischaemic Stroke")

x$Predictor <- str_replace(x$Predictor, "gdf15", "GDF15")
x$Predictor <- str_replace(x$Predictor, "nt.probnp", "Nt-proBNP")

# x$Outcome2 <- x$Outcome
# x$TraitVar <- paste0(x$Name)
# x$TraitVar = factor(x$TraitVar, levels=unique(x$TraitVar[rev(order(x$Hazard.Ratio.x))]))

My_Theme = theme(
  panel.border = element_rect(colour="black",size=1, fill = NA),
  axis.title.x = element_text(size = 20), # controls HR label size 
  axis.text.x = element_text(size = 20),
  axis.text.y = element_text(size = 20),
  axis.title.y = element_text(size = 20),
  strip.text = element_text(size = 20, face = "bold"),
  legend.text=element_text(size=20),
  legend.title=element_text(size=20, face = "bold"), legend.position = "none",
  axis.title=element_text(size=20))


# collist1 <- c("yellow", "blue", "turquoise", "green", "purple", "pink", "grey", "orange")
# collist2 <- c("yellow", "blue", "turquoise", "green", "purple", "pink", "grey", "orange")

collist1 <- c("red", "black", "red", "red", "red", "red", "black", "red")
collist2 <- c("red", "black", "black", "red", "red", "red", "red", "red")

x$Predictor <- str_replace(x$Predictor, "Nt-proBNP", "NT-proBNP")

pdf("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/01_Cox_20k_results/20k_plot.pdf", width = 18, height = 5)
ggplot(x,aes(y=Hazard.Ratio, x=Predictor)) + 
  geom_point(size = 4.5, colour = collist1)+
  geom_errorbar(aes(ymin = LCI, ymax = UCI),
                position = position_dodge(0.5), width = 0.1,
                colour = collist2)+
ylab("Hazard Ratio (95% CI)")+ xlab ("") + theme_classic() +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme
  dev.off()
