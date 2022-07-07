library(atc)
library(ggplot2)

atcdf <- atcs(atc, filter=LevelFilter(5))
total_fraction_atc <- substr(atcdf$key, start = 1, stop = 1)

json_opentargets_parsed <- readRDS("/path/to/final_opentargets_parsed.rds")

drugable_gene_modules <- list()

for(i in 1:length(members_list)){
  df <- json_opentargets_parsed[which(json_opentargets_parsed$dat.target.gene_info.symbol %in% members_list[[i]]),]
  drugable_gene_modules[[i]] <- df
}

drug_statistics <- data.frame(module_dim=sapply(members_list, length), num_genes=sapply(drugable_gene_modules, function(x) length(unique(x$dat.target.gene_info.symbol))), 
                              num_drugs=sapply(drugable_gene_modules, function(x) length(unique(x$dat.drug.id))),
                              module=paste0("mod_", rownames(drug_statistics)))
drug_statistics$module <- NULL

df <- reshape::melt(t(drug_statistics))
colnames(df) <- c("feature", "module", "value")

ggplot(df, aes(fill=feature, y=value, x=reorder(module))) + 
  geom_bar(position="dodge", stat="identity") + ggtitle("") + theme_minimal()
#num_associations=sapply(drugable_gene_modules, nrow)

plot_atc_values <- function(drug_list, mode="fraction", title=""){
  atcdf <- atc::atcs(atc, filter=LevelFilter(5))
  total_fraction_atc <- substr(atcdf$key, start = 1, stop = 3)
  
  atc_BRCA <- data.frame(drug_list)
  atc_BRCA$atc <- NA
  atc_BRCA$atc <- atcdf$key[match(tolower(atc_BRCA$drug_list), atcdf$name)]
  atc_BRCA$atc <- substr(atc_BRCA$atc, start = 1, stop = 3)
  
  fraction_atc_BRCA <- table(atc_BRCA$atc)
  fraction_total <- table(total_fraction_atc)
  fraction_total <- fraction_total[which(names(fraction_total) %in% names(fraction_atc_BRCA))]
  gdata_BRCA <- data.frame(fraction_atc_BRCA, fraction_total)
  rownames(gdata_BRCA) <- names(fraction_total)
  gdata_BRCA$Var1 <- gdata_BRCA$total_fraction_atc <- NULL
  gdata_BRCA$value <- (gdata_BRCA$Freq/gdata_BRCA$Freq.1)*100
  gdata_BRCA$element <- rownames(gdata_BRCA)
  if(mode=="fraction"){
    p <- ggplot(gdata_BRCA, aes(fill=element, y=value, x=reorder(element, -value))) + 
      geom_bar(position="dodge", stat="identity") + ggtitle(title) + theme_minimal()
    
    p
    
  }else if(mode=="absolute"){
    p <- ggplot(gdata_BRCA, aes(fill=element, y=Freq, x=reorder(element, -value))) + 
      geom_bar(position="dodge", stat="identity")+ ggtitle(title) + theme_minimal()
    
    p
  }
  
  
}



mod1_atc <- plot_atc_values(tolower(drugable_gene_modules[[1]]$dat.drug.molecule_name))
mod2_atc <- plot_atc_values(tolower(drugable_gene_modules[[2]]$dat.drug.molecule_name))
mod3_atc <- plot_atc_values(tolower(drugable_gene_modules[[3]]$dat.drug.molecule_name))
mod4_atc <- plot_atc_values(tolower(drugable_gene_modules[[4]]$dat.drug.molecule_name))
mod5_atc <- plot_atc_values(tolower(drugable_gene_modules[[5]]$dat.drug.molecule_name))
mod6_atc <- plot_atc_values(tolower(drugable_gene_modules[[6]]$dat.drug.molecule_name))
mod7_atc <- plot_atc_values(tolower(drugable_gene_modules[[7]]$dat.drug.molecule_name))
mod8_atc <- plot_atc_values(tolower(drugable_gene_modules[[8]]$dat.drug.molecule_name))
mod9_atc <- plot_atc_values(tolower(drugable_gene_modules[[9]]$dat.drug.molecule_name))
mod11_atc <- plot_atc_values(tolower(drugable_gene_modules[[11]]$dat.drug.molecule_name))


###### DrugBank ######

setwd("/path/to/drugbank/")
vocab <- read_csv("drugbank vocabulary.csv", col_names = TRUE)

all.files <- list.files(".", pattern = "all.csv", recursive = TRUE)

drugbank_list <- list()

for(i in 1:length(all.files)){
  df1 <- read_csv(all.files[i], col_names = TRUE)
  drugbank_list[[i]] <- df1
  
}

all_drugbank <- do.call("rbind", drugbank_list)

all_drugbank <- all_drugbank %>% separate_rows(`Drug IDs`)

all_drugbank$drug_name <- NA
all_drugbank$drug_name <- vocab$`Common name`[match(all_drugbank$`Drug IDs`, vocab$`DrugBank ID`)]

drugable_gene_modules <- list()

for(i in 1:length(members_list)){
  df <- all_drugbank[which(all_drugbank$`Gene Name` %in% members_list[[i]]),]
  drugable_gene_modules[[i]] <- df
}

drug_statistics <- data.frame(module_dim=sapply(members_list, length), num_genes=sapply(drugable_gene_modules, function(x) length(unique(x$`Gene Name`))), 
                              num_drugs=sapply(drugable_gene_modules, function(x) length(unique(x$`Drug IDs`))),
                              module=names(members_list))
drug_statistics$module <- NULL

drug_statistics$freq_druggable <- drug_statistics$num_genes/drug_statistics$module_dim

df <- reshape::melt(t(drug_statistics))
colnames(df) <- c("feature", "module", "value")

ggplot(df, aes(fill=feature, y=value, x=reorder(module))) + 
  geom_bar(position="dodge", stat="identity") + ggtitle("") + theme_minimal()
#num_associations=sapply(drugable_gene_modules, nrow)



mod1_atc <- plot_atc_values(tolower(drugable_gene_modules[[1]]$drug_name), mode = "fraction")
mod2_atc <- plot_atc_values(tolower(drugable_gene_modules[[2]]$drug_name), mode = "fraction")
mod3_atc <- plot_atc_values(tolower(drugable_gene_modules[[3]]$drug_name), mode = "fraction")
mod4_atc <- plot_atc_values(tolower(drugable_gene_modules[[4]]$drug_name), mode = "fraction")
mod5_atc <- plot_atc_values(tolower(drugable_gene_modules[[5]]$drug_name), mode = "fraction")
mod6_atc <- plot_atc_values(tolower(drugable_gene_modules[[6]]$drug_name), mode = "fraction")
mod7_atc <- plot_atc_values(tolower(drugable_gene_modules[[7]]$drug_name), mode = "fraction")
mod8_atc <- plot_atc_values(tolower(drugable_gene_modules[[8]]$drug_name), mode = "fraction")
mod9_atc <- plot_atc_values(tolower(drugable_gene_modules[[9]]$drug_name), mode = "fraction")
mod11_atc <- plot_atc_values(tolower(drugable_gene_modules[[11]]$drug_name), mode = "fraction")

for (i in 1:length(drugable_gene_modules)){
  drugable_gene_modules[[i]]$module <- paste0("mod_", i)
}

total_drugability_df <- do.call("rbind", drugable_gene_modules)

drugab_freq <- as.data.frame(table(total_drugability_df$drug_name))
drugab_freq$num_mod <- NA
drugab_freq$Var1 <- as.character(drugab_freq$Var1)

for(i in 1:length(drugab_freq$Var1)){
  
  df_drug <- total_drugability_df[which(total_drugability_df$drug_name %in% drugab_freq$Var1[[i]]),]
  drugab_freq$num_mod[i] <- length(table(df_drug$module))
  
}

drugab_freq_1_mod <- drugab_freq[which(drugab_freq$num_mod==1),]
drugab_freq_1_mod$mod <- NA
drugab_freq_1_mod$mod <- total_drugability_df$module[match(drugab_freq_1_mod$Var1, total_drugability_df$drug_name)]

table(drugab_freq_1_mod$mod)


drug_1_mod_freq <- data.frame(module_dim=drug_statistics$module_dim, drugs=NA)
rownames(drug_1_mod_freq) <- rownames(drug_statistics)
drug_1_mod_freq$drugs <- table(drugab_freq_1_mod$mod)[match(rownames(drug_1_mod_freq), names(table(drugab_freq_1_mod$mod)))]
drug_1_mod_freq$freq <- drug_1_mod_freq$drugs/drug_1_mod_freq$module_dim

mod1_atc_1 <- plot_atc_values(tolower(drugab_freq_1_mod$Var1[which(drugab_freq_1_mod$mod=="mod_1")]), mode = "fraction")
mod2_atc_1 <- plot_atc_values(tolower(drugab_freq_1_mod$Var1[which(drugab_freq_1_mod$mod=="mod_2")]), mode = "fraction")
mod3_atc_1 <- plot_atc_values(tolower(drugab_freq_1_mod$Var1[which(drugab_freq_1_mod$mod=="mod_3")]), mode = "fraction")
mod4_atc_1 <- plot_atc_values(tolower(drugab_freq_1_mod$Var1[which(drugab_freq_1_mod$mod=="mod_4")]), mode = "fraction")
mod5_atc_1 <- plot_atc_values(tolower(drugab_freq_1_mod$Var1[which(drugab_freq_1_mod$mod=="mod_5")]), mode = "fraction")
mod6_atc_1 <- plot_atc_values(tolower(drugab_freq_1_mod$Var1[which(drugab_freq_1_mod$mod=="mod_6")]), mode = "fraction")
mod7_atc_1 <- plot_atc_values(tolower(drugab_freq_1_mod$Var1[which(drugab_freq_1_mod$mod=="mod_7")]), mode = "fraction")
mod8_atc_1 <- plot_atc_values(tolower(drugab_freq_1_mod$Var1[which(drugab_freq_1_mod$mod=="mod_8")]), mode = "fraction")
mod9_atc_1 <- plot_atc_values(tolower(drugab_freq_1_mod$Var1[which(drugab_freq_1_mod$mod=="mod_9")]), mode = "fraction")
mod11_atc_1 <- plot_atc_values(tolower(drugab_freq_1_mod$Var1[which(drugab_freq_1_mod$mod=="mod_11")]), mode = "fraction")


total_drugability_df$atc <- NA
total_drugability_df$atc <- atcdf$key[match(tolower(total_drugability_df$drug_name), atcdf$name)]
tmp <- total_drugability_df[which(total_drugability_df$module=="mod_8"),]
tmp$`GeneCard ID` <- NULL
tmp <- na.omit(tmp)
tmp$rank <- NA
tmp$rank <- named_ranklist[match(tmp$`Gene Name`, names(named_ranklist))]



atcdf <- atc::atcs(atc, filter=LevelFilter(5))
total_fraction_atc <- substr(atcdf$key, start = 1, stop = 3)

atc_BRCA <- data.frame(tolower(drugab_freq_1_mod$Var1[which(drugab_freq_1_mod$mod %in% c("mod_2", "mod_8", "mod_11", "mod_4"))]))
colnames(atc_BRCA) <- "drug_list"
atc_BRCA$atc <- NA
atc_BRCA$atc <- atcdf$key[match(tolower(atc_BRCA$drug_list), atcdf$name)]
atc_BRCA$atc <- substr(atc_BRCA$atc, start = 1, stop = 3)

fraction_atc_BRCA <- table(atc_BRCA$atc)
fraction_total <- table(total_fraction_atc)
fraction_total <- fraction_total[which(names(fraction_total) %in% names(fraction_atc_BRCA))]
gdata_BRCA <- data.frame(fraction_atc_BRCA, fraction_total)
rownames(gdata_BRCA) <- names(fraction_total)
gdata_BRCA$Var1 <- gdata_BRCA$total_fraction_atc <- NULL
gdata_BRCA$value <- (gdata_BRCA$Freq/gdata_BRCA$Freq.1)*100
gdata_BRCA$element <- rownames(gdata_BRCA)
gdata_BRCA$letter <- substr(gdata_BRCA$element, start = 1, stop = 1)
p <- ggplot(gdata_BRCA, aes(fill=letter, y=value, x=reorder(element, -value))) + 
  geom_bar(position="dodge", stat="identity") + theme_minimal() + scale_fill_manual(values = c("A" = "aquamarine2", "B" = "blue", "C"="coral", "D" = "blueviolet", "G"="brown2", "J"="darkgrey", "L" = "cadetblue2", "N"="darkgreen", "P" = "darksalmon", "R"="darkorange", "S"="cornflowerblue", "V"= "darkolivegreen2", "M"="cornsilk3")) +
 scale_x_discrete(name ="ATC category") + theme(axis.text.x = element_text(size=14),
                                                  axis.text.y = element_text(size=12))

p
