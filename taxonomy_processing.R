taxonomy_16S_step1 <-
  read_delim("~/ktaxonomy_manual.txt", col_names = F) %>%
  select(1, 3, 5, 9) %>%
  `colnames<-`(c("tax_id", "parent", "type", "name")) %>%
  column_to_rownames("tax_id")
for (x in 1:nrow(taxonomy_16S_step1)) {
  if (x == 1) {
    taxonomy_16S_step1[x, "taxonomy"] <- "root"
  } else {
    taxonomy_16S_step1[x, "taxonomy"] <- paste0(taxonomy_16S_step1[as.character(taxonomy_16S_step1[x, "parent"]), "taxonomy"], ";", taxonomy_16S_step1[x, "name"])
  }
}
taxonomy_16S <-
  taxonomy_16S_step1 %>%
  rownames_to_column("tax_id") %>%
  filter(tax_id %in% unique(c(s16_kraken$tax_id, kraken_raw$tax_id, rna$tax_id))) %>%
  mutate(type = str_replace(type, "R", "root")) %>%
  mutate(type = str_replace(type, "D", "superkingdom")) %>%
  mutate(type = str_replace(type, "P", "phylum")) %>%
  mutate(type = str_replace(type, "C", "class")) %>%
  mutate(type = str_replace(type, "O", "order")) %>%
  mutate(type = str_replace(type, "F", "family")) %>%
  mutate(type = str_replace(type, "G", "genus")) %>%
  mutate(type = str_replace(type, "S", "species")) %>%
  filter(tax_id %in% unique(c(s16_kraken$tax_id, kraken_raw$tax_id, rna$tax_id))) %>%
  select(-parent) %>%
  as_tibble()
rm(taxonomy_16S_step1)


#### RNAseq taxonomy
taxonomy_rna_step1 <-
  read_delim("~/ktaxonomy_rna.tsv", col_names = F) %>%
  select(1, 3, 5, 9) %>%
  `colnames<-`(c("tax_id", "parent", "type", "name")) %>%
  column_to_rownames("tax_id")
for (x in 1:nrow(taxonomy_rna_step1)) {
  if (x == 1) {
    taxonomy_rna_step1[x, "taxonomy"] <- "root"
  } else {
    taxonomy_rna_step1[x, "taxonomy"] <- paste0(taxonomy_rna_step1[as.character(taxonomy_rna_step1[x, "parent"]), "taxonomy"], ";", taxonomy_rna_step1[x, "name"])
  }
}
taxonomy_RNA <-
  taxonomy_rna_step1 %>%
  rownames_to_column("tax_id") %>%
  filter(tax_id %in% unique(c(s16_kraken$tax_id, kraken_raw$tax_id, rna$tax_id))) %>%
  mutate(type = str_replace(type, "R", "root")) %>%
  mutate(type = str_replace(type, "D", "superkingdom")) %>%
  mutate(type = str_replace(type, "P", "phylum")) %>%
  mutate(type = str_replace(type, "C", "class")) %>%
  mutate(type = str_replace(type, "O", "order")) %>%
  mutate(type = str_replace(type, "F", "family")) %>%
  mutate(type = str_replace(type, "G", "genus")) %>%
  mutate(type = str_replace(type, "S", "species")) %>%
  filter(tax_id %in% unique(c(s16_kraken$tax_id, kraken_raw$tax_id, rna$tax_id))) %>%
  select(-parent) %>%
  as_tibble()
rm(taxonomy_rna_step1)

#### RNAseq taxonomy
taxonomy_wgs_step1 <-
  read_delim("~/ktaxonomy_rna.tsv", col_names = F) %>%
  select(1, 3, 5, 9) %>%
  `colnames<-`(c("tax_id", "parent", "type", "name")) %>%
  column_to_rownames("tax_id")
for (x in 1:nrow(taxonomy_wgs_step1)) {
  if (x == 1) {
    taxonomy_wgs_step1[x, "taxonomy"] <- "root"
  } else {
    taxonomy_wgs_step1[x, "taxonomy"] <- paste0(taxonomy_wgs_step1[as.character(taxonomy_wgs_step1[x, "parent"]), "taxonomy"], ";", taxonomy_wgs_step1[x, "name"])
  }
}
taxonomy_WGS <-
  taxonomy_wgs_step1 %>%
  rownames_to_column("tax_id") %>%
  filter(tax_id %in% unique(c(s16_kraken$tax_id, kraken_raw$tax_id, rna$tax_id))) %>%
  mutate(type = str_replace(type, "R", "root")) %>%
  mutate(type = str_replace(type, "D", "superkingdom")) %>%
  mutate(type = str_replace(type, "P", "phylum")) %>%
  mutate(type = str_replace(type, "C", "class")) %>%
  mutate(type = str_replace(type, "O", "order")) %>%
  mutate(type = str_replace(type, "F", "family")) %>%
  mutate(type = str_replace(type, "G", "genus")) %>%
  mutate(type = str_replace(type, "S", "species")) %>%
  filter(tax_id %in% unique(c(s16_kraken$tax_id, kraken_raw$tax_id, rna$tax_id))) %>%
  select(-parent) %>%
  as_tibble()
rm(taxonomy_wgs_step1)

kraken_taxonomy<-
  taxonomy_RNA %>%
  rbind(taxonomy_WGS %>% filter(!tax_id %in% taxonomy_RNA$tax_id))
unified_taxonomy <-
  kraken_taxonomy%>%
  rbind(taxonomy_16S %>% filter(!tax_id %in% taxonomy_RNA$tax_id))
kraken_microbeTaxonomy <-
  unified_taxonomy %>%
  filter(str_detect(taxonomy, "Bacteria"))

all_bact_genera <-
  kraken_microbeTaxonomy %>%
  filter(type == "genus") %>%
  pull(tax_id)
all_bact_species <-
  kraken_microbeTaxonomy %>%
  filter(type == "species") %>%
  pull(tax_id)
all_bact_phyla <-
  kraken_microbeTaxonomy %>%
  filter(type == "phylum") %>%
  pull(tax_id)
all_bact_family <-
  kraken_microbeTaxonomy %>%
  filter(type == "family") %>%
  pull(tax_id)

# read taxonomy and limit to only the bacteria in the dataset
relevant_taxa<-
  read_delim('extra_data/ncbi_taxonomy.txt')%>%
  filter(tax_id %in% union(rna$tax_id,
                           union(s16_kraken$tax_id,
                                 kraken_raw$tax_id)))%>%
  mutate(tax_id=as.character(tax_id))
tax_table<-build_taxTable(relevant_taxa)