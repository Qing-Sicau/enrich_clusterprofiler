library(tidyr)
library(dplyr)

up_gene2k <- read.csv("upregulated_gene2k", header = F, sep = "\t")
colnames(up_gene2k) <- c("gene_id","k_number")
head(up_gene2k)

k2komap <- read.csv("complete_knumber_2_pathway.txt", header = T, sep = "\t")
head(k2komap)

up_gene2komap <- up_gene2k %>%
  dplyr::left_join(k2komap, by = c("k_number" = "k_number")) %>%
  dplyr::select(gene_id, komap) %>%
  na.omit()
head(up_gene2komap)

plant_sp_komap2desc <- read.csv("plant_specific_kegg_map.txt", header = T, sep = "\t")
head(plant_sp_komap2desc)

up_gene2komap_desc <- up_gene2komap %>%
  dplyr::left_join(plant_sp_komap2desc, by = c("komap" = "Pathway")) %>%
  dplyr::select(gene_id, Name) %>%
  na.omit()

up_gene2kegg <- up_gene2komap_desc %>% 
  group_by(Name) %>% 
  summarize(number = n_distinct(gene_id)) %>% 
  arrange(desc(number)) %>%
  top_n(20)

write.csv(up_gene2kegg, "up_gene2keegg.csv", quote = F)

down_gene2k <- read.csv("down_gene2k", header = F, sep = "\t")
colnames(down_gene2k) <- c("gene_id","k_number")
head(down_gene2k)

k2komap <- read.csv("complete_knumber_2_pathway.txt", header = T, sep = "\t")
head(k2komap)

down_gene2komap <- down_gene2k %>%
  dplyr::left_join(k2komap, by = c("k_number" = "k_number")) %>%
  dplyr::select(gene_id, komap) %>%
  na.omit()
head(down_gene2komap)

plant_sp_komap2desc <- read.csv("plant_specific_kegg_map.txt", header = T, sep = "\t")
head(plant_sp_komap2desc)

down_gene2komap_desc <- down_gene2komap %>%
  dplyr::left_join(plant_sp_komap2desc, by = c("komap" = "Pathway")) %>%
  dplyr::select(gene_id, Name) %>%
  na.omit()

down_gene2kegg <- down_gene2komap_desc %>% 
  group_by(Name) %>% 
  summarize(number = n_distinct(gene_id)) %>% 
  arrange(desc(number)) %>%
  top_n(20)

write.csv(down_gene2kegg, "down_gene2keegg.csv", quote = F)