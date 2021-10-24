# For local db construction and enrichment analysis
# updated 2021-10-24;
# supnovel@sicau.edu.cn
# the json section canbe adjusted to accommodate for other rel.

# setwd("E:/4-RNA-seq/fra_ann_enrich")

# 1- load packages needed
library("dplyr")
library("stringr")
library("tidyr")
library("clusterProfiler")
# BiocManager::install("clusterProfiler")

# 2- load data
myGeneData <- read.csv("gene2gi.tsv", header = T, sep = "\t", stringsAsFactors = F)
head(myGeneData)
dim(myGeneData)

## 2.1 gene id and description
fa_g2id <- myGeneData %>%
  dplyr::select(GID = seq_id, NT_seq = nt_seq_id, REF = hit_id, DESCRIPTION = desc)
head(fa_g2id)
dim(fa_g2id)

## 2.2 gene and go data
## TOA pipeline results were piped and selected using sed and awk.
myGene2GoData <- read.csv("gene2go.ok.tsv", header = T, sep = "\t", stringsAsFactors = F)
dim(myGene2GoData)
fa_g2go <- myGene2GoData %>%
  dplyr::select(GID = gene_id, GO = go_id) %>%
  dplyr::mutate(EVIDENCE = "IEA")
head(fa_g2go)
dim(fa_g2go)

## 2.3 gene and kegg data
## KAAS annotation results, sort by k_number, remove all genes that no pathways were assigned.

myGene2KData <- read.csv("gene2kegg.tsv", header = T, sep = "\t", stringsAsFactors = F)
dim(myGene2KData)
fa_g2k <- myGene2KData %>%
  dplyr::select(GID = gene_id, KNUMBER = K_number)
dim(fa_g2k)

## 2.4 check the header
head(fa_g2id,1)
head(fa_g2go,1)
head(fa_g2k,1)

## 2.5 remove all duplicates
gene_info <- fa_g2id %>%
  na.omit()
gene2go <- fa_g2go %>%
  na.omit()
gene2k <- fa_g2k %>%
  na.omit()

## 2.6 the correspondence of k to ko(map), ko(map) to description
## this can be obtained in various way

## 2.6.1 automatically obtain through ko00001 of the kegg site
## https://www.genome.jp/kegg-bin/get_htext?ko00001
# library(jsonlite)
# library(purrr)
# library(RCurl)
# run only once when needed
# if(T){
#   update_kegg <- function(json = "ko00001.json"){
#     komap2description <- tibble(komap = character(), description = character())
#     k2komap <- tibble(k_number = character(), komap = character())
#     kegg <- fromJSON(json)
#     for (a in seq_along(kegg[["children"]][["children"]])){
#       A <- kegg[["children"]][["name"]][[a]]
#       for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])){
#         B <- kegg[["children"]][["children"]][[a]][["name"]][[b]]
#         for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])){
#           komap_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
#           komap_id <- str_match(komap_info, "ko[0-9]{5}")[1]
#           komap_description <- str_replace(komap_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
#           komap2description <- rbind(komap2description, tibble(komap = komap_id, description = komap_description))
#           kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
#           kos <- str_match(kos_info, "K[0-9]*")[,1]
#           k2komap <- rbind(k2komap, tibble(k_number = kos, komap = rep(komap_id, length(kos))))
#         }
#       }
#     }
#     save(komap2description, k2komap, file = "kegg_info.RData")
#   }
#   update_kegg(json = "ko00001.json")
# }

# each time of run, direct start here
# load(file = "kegg_info.RData")

## 2.6.2 directly download from the site
# http://rest.kegg.jp/list/pathway/

komap2description_data <- read.csv("complete_keggmap2description.txt", header = F, sep = "\t")
komap2description <- data.frame(lapply(komap2description_data, function(x){str_replace(x, "path:map", "ko")}))
colnames(komap2description) <- c("komap", "description")

# http://rest.kegg.jp/link/pathway/ko
# https://www.kegg.jp/kegg/rest/keggapi.html
# more link informatin could be found here, for k2ko, ko2description, k2module, ko2module etc.

k2komap_data <- read.csv("complete_knumber2keggmap.txt", header = F, sep = "\t")
k2komap <- k2komap_data %>%
  dplyr::filter(grepl("path:ko", V2)) %>%
  tidyr::separate(V1, into = c("x", "k_number"), sep = ":") %>%
  tidyr::separate(V2, into = c("y","komap"), sep = ":") %>%
  dplyr::select(k_number, komap) %>%
  unique()

## 2.6.3 gene to komap and descripton
head(gene2k, 5)
head(k2komap, 5)

gene2komap <- gene2k %>%
  dplyr::left_join(k2komap, by = c("KNUMBER" = "k_number")) %>%
  dplyr::select(GID, komap) %>%
  na.omit()
head(gene2komap)

## 2.6.4 plant specific komaps
komap_plant <- read.csv("plant_specific_kegg_map.txt", header = T, sep = "\t")

## 2.6.5 only kept the annotated pathway of plant in the gene list
head(gene2komap)
gene2komap_plant <- gene2komap %>%
  dplyr::inner_join(komap_plant, by = c("komap" = "Pathway"))

# 3 construct the local database
#BiocManager::install('AnnotationForge')
library(AnnotationForge)

## search for the Taxonomy，for example fragaria
## https://www.ncbi.nlm.nih.gov/taxonomy/?term=

tax_id = "3747"
genus = "Fragaria"
species = "ananassa"
nrow(gene_info)
gene_info<-gene_info[!duplicated(gene_info),]
nrow(gene_info)
nrow(gene2go)
gene2go<-gene2go[!duplicated(gene2go),]
nrow(gene2go)

nrow(gene2k)
gene2k<-gene2k[!duplicated(gene2k),]
nrow(gene2k)

nrow(gene2komap_plant)
gene2komap_plant <- gene2komap_plant[!duplicated(gene2komap_plant),]
nrow(gene2komap_plant)

colnames(gene_info)
colnames(gene2go)
colnames(gene2k)
colnames(gene2komap_plant)

#
# gene2go must include three columns, and must be GID, GOID and EVIDENCE respectively
# 

makeOrgPackage(gene_info=gene_info,
               go=gene2go,
               ko=gene2k,
               maintainer = "supnovel@sicau.edu.cn",
               author = "Q.Chen",
               pathway=gene2komap_plant,
               version="1.02a",
               outputDir = ".",
               tax_id=tax_id,
               genus=genus,
               species=species,
               goTable="go")


# 4 We can now use this db to do enrichment analysis

library(AnnotationHub)
library(AnnotationDbi)
library(clusterProfiler)
library("dplyr")
library("stringr")
library("tidyr")

fraaDb <- loadDb("org.Fananassa.eg.sqlite")
k2pathway <- read.csv("complete_knumber_2_pathway.txt", header = T, sep = "\t")
pathway2name <- read.csv("plant_specific_kegg_map.txt", header = T, sep = "\t")

# ?AnnotationDb
# columns(x)
# keytypes(x)
# keys(x, keytype, ...)
# select(x, keys, columns, keytype, ...)
# mapIds(x, keys, column, keytype, ..., multiVals)
# saveDb(x, file)
# loadDb(file, packageName=NA)

columns(fraaDb)
keytypes(fraaDb)
head(keys(fraaDb, keytype="GID"),10)

myDEGs <- read.csv("dgeResultExactTest.csv", header = T, sep = ",", stringsAsFactors = F)
head(myDEGs,5)

myDEGs_list <- myDEGs %>%
  dplyr::filter(grepl("\\d", logFC, perl = T)) %>%
  dplyr::distinct(transcript_id, .keep_all = T) %>%
  dplyr::filter(abs(logFC) >= 0.5 & FDR <= 0.05)

## for GO and KEGG enrich
head(myDEGs_list)
gene_list <- myDEGs_list$transcript_id

## the term pathway2gene can not be changed here
pathway2gene <- AnnotationDbi::select(fraaDb, 
                                      keys = keys(fraaDb), 
                                      columns = c("komap", "GID")) %>%
              na.omit() %>% 
              dplyr::select(komap, GID)

## the term patway2name can not be changed here
pathway2name <- AnnotationDbi::select(fraaDb,
                                           keys = keys(fraaDb),
                                           columns = c("komap", "Name")) %>%
                    na.omit() %>%
                    dplyr::select(komap, Name)

## 4.1 kegg pathway enrichment 
enrichedKEEGpathway <- enricher(gene = gene_list, 
                                pvalueCutoff = 0.1, 
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.5, 
                                TERM2GENE = pathway2gene,
                                TERM2NAME = pathway2name)

barplot(enrichedKEEGpathway, showCategory=20,color="pvalue", font.size=10)

## 4.2 go term enrichment
my_enrichGO <- enrichGO(gene = gene_list,
                        OrgDb = fraaDb,
                        keyType = "GID",
                        ont = "MF",
                        pvalueCutoff = 0.1, 
                        pAdjustMethod = "BH",
                        minGSSize = 10, readable = F)

barplot(my_enrichGO, showCategory=20, x = "GeneRatio")

# 5 GSEA analysis
myDEGs_gsea <- read.csv("dgeResultExactTest.csv", header = T, sep = ",", stringsAsFactors = F)
head(myDEGs_gsea,5)

myDEGs_gsea_list <- myDEGs_gsea %>%
  dplyr::filter(grepl("\\d", logFC, perl = T)) %>%
  dplyr::distinct(transcript_id, .keep_all = T)

## for GO and KEGG GSEA
head(myDEGs_gsea_list)
gene_gsea_list <- myDEGs_gsea_list$logFC
names(gene_gsea_list) <- myDEGs_gsea_list$transcript_id
gene_gsea_list <- sort(gene_gsea_list, decreasing = T)
head(gene_gsea_list)

gse_CC <- gseGO(geneList = gene_gsea_list, 
                ont = "mf", 
                OrgDb = fraaDb, 
                keyType = "GID", 
                nPerm = 1000, 
                minGSSize = 10, 
                pvalueCutoff = 0.1, 
                pAdjustMethod = "BH")

head(summary(gse_CC))
gseaplot(gse_CC, geneSetID = "GO:0004857")

gene_map <- bitr(myDEGs_gsea_list$transcript_id, 
                 fromType = "GID", 
                 toType = "komap", 
                 OrgDb = fraaDb)

gene_map <- gene_map %>%
  dplyr::distinct(GID, .keep_all=T) %>%
  na.omit()

head(gene_map)
term2gene <- gene_map[,c(2,1)]

gse_kegg <- GSEA(geneList = gene_gsea_list, 
                 #nPerm = 3000, 
                 pvalueCutoff = 0.05, 
                 pAdjustMethod = "fdr", 
                 TERM2GENE = term2gene, 
                 minGSSize = 10,
                 TERM2NAME = pathway2name)

GSEA_results_kegg <- gse_kegg@result
write.csv(GSEA_results_kegg, "gsea_results_kegg.csv")

## Plotting
library(enrichplot)
library(ggplot2)
library("scales")
library(RColorBrewer)
clusterProfiler::dotplot(gse_kegg, 
                         showCategory = 10, 
                         split = ".sign")+ facet_grid(~.sign)

enrichplot::gseaplot2(x = gse_kegg, geneSetID = c("ko04141","ko00250"), pvalue_table = T)

sort_gsea_kegg <- gse_kegg[order(GSEA_results_kegg$enrichmentScore, decreasing =T),]
head(sort_gsea_kegg)

gseaplot2(gse_kegg,
          row.names(sort_gsea_kegg)[1:3],
          color = colorspace::rainbow_hcl(1), 
          pvalue_table = T)

## plot using ggplot2
ggplot(sort_gsea_kegg[c(1:10),],aes(x = NES,
                                    y = Description,
                                    color = p.adjust,
                                    size = setSize))+
  geom_point()+
  scale_size_continuous(breaks = c(10,20,40,80,160),
                        name = "No. of\nsignificant DEGs")+
  scale_color_gradient2(low = "black", mid = "blue", high = "red",
                        limits=c(0,0.1),
                        breaks = seq(0.02, 0.08, 0.02),
                        labels = seq(0.02, 0.08, 0.02))+
  xlim(c(-3,3))+
  geom_vline(xintercept=0, linetype = 2, color = "gray")+
  xlab("Normalized enrichment score")+
  ylab("")+
  theme_bw()+
  guides(color = guide_colorbar(reverse = T))

myDEGs <- read.csv("dgeResultExactTest.csv", header = T, sep = ",", stringsAsFactors = F)
head(myDEGs,5)

myDEGs_list <- myDEGs %>%
  dplyr::filter(grepl("\\d", logFC, perl = T)) %>%
  dplyr::distinct(transcript_id, .keep_all = T) %>%
  dplyr::filter(abs(logFC) >= 0.5 & FDR <= 0.05)

## for GO and KEGG enrich
head(myDEGs_list)

dge2k <- read.csv("DEG-Ko.txt", header = F, sep = "\t" , fill = T)
colnames(dge2k) <- c("transcript_id", "k_number")

dge2knumber <- dge2k %>%
  dplyr::select(transcript_id, k_number) %>%
  dplyr::filter(k_number != "")

knumber2module <- read.csv("complete_knumber_2_modules.txt", header = F, sep = "\t", stringsAsFactors = F)
colnames(knumber2module) <- c("k_number", "module_id")

knumber2module_m <- data.frame(lapply(knumber2module, function(x){str_replace(x, "ko:", "")}))
knumber2module_ok <- data.frame(lapply(knumber2module_m, function(x){str_replace(x, "md:", "")})) 

module2description <- read.csv("complete_module_2_description.txt", header = F, sep = "\t", strip.white = F)
colnames(module2descripton) <- c("module_id", "description")
module2description <- data.frame(lapply(module2descripton, function(x){str_remove(x, "md:")}))

head(dge2knumber)
head(knumber2module_ok)

gene2module <- dge2knumber %>%
  dplyr::left_join(knumber2module_ok, by = c("k_number" = "k_number")) %>%
  dplyr::select(transcript_id, module_id) %>%
  na.omit() %>%
  dplyr::left_join(module2description, by = c("module_id" = "module_id"))
write.table(gene2module, "dge_kegg_modules.txt", row.names = F, quote = F, sep = "\t")

dge_kegg_module_stats <- gene2module %>%
  dplyr::group_by(module_id) %>%
  dplyr::summarise(number = n_distinct(transcript_id)) %>%
  dplyr::arrange(desc(number)) %>%
  dplyr::left_join(module2description, by = c("module_id" = "module_id"))
write.table(dge_kegg_module_stats, "dge_kegg_module_statistics.txt", sep = "\t", row.names = F, quote = F)