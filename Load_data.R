source("microbiome_functions.R")

#############################################################################################
############################### READ IN DATA #################################################
#############################################################################################
load("~/Documents/GitHub/Sherlock-Lung-Data-Analysis/Data/clinical_data.RData")
load("~/Documents/GitHub/Sherlock-Lung-Data-Analysis/Data/BBsolution_final3_short.RData")
load("~/Documents/GitHub/Sherlock-Lung-Data-Analysis/Data/sherlock_data_all.RData")
load("~/Downloads/SBScolor.RData")
load("extra_data/sherlock_sigantures.RData")
SBScolor<-unname(SBScolor)
library(RColorBrewer)
brewer.pal(3,name = "Accent")
s16_color<-brewer.pal(3,name = "Accent")[1]
RNA_color<-brewer.pal(3,name = "Accent")[2]
WGS_color<-brewer.pal(3,name = "Accent")[3]

# some annotations to load
broad_big_list <-
  read_delim("~/Downloads/wgs_broad_biglist.csv", delim = ",")
harmonization <-
  read_delim("~/Downloads/Sherlock_Harmonization_Final_112221.txt") %>%
  mutate(Sherlock_PID = if_else(is.na(Sherlock_PID), PATIENT_ID, Sherlock_PID)) %>%
  mutate(Sherlock_PID = str_replace(Sherlock_PID, "NSLC-", ""))
sherlock_tp53_status<-
  sherlock_data_full%>%
  filter(Gene=="TP53",Type=="Mutation_Driver")%>%
  select(Subject,Alteration)
  

# WGS annotations; 
annotations <-
  bind_rows(
    wgs_groups_info %>% select(Subject, Barcode = Tumor_Barcode, Study, Assigned_Population:SP_Group) %>% mutate(Sample_Source = "Tumor"),
    wgs_groups_info %>% select(Subject, Barcode = Normal_Barcode, Study, Assigned_Population:SP_Group, Sample_Source = Normal_Material)
  ) %>%
  select(-SBS4)
annotations_withClinical <-
  full_join(annotations, clinical_data) %>%
  mutate(Stage_simple = str_replace(Stage, "[AB].*$", ""))
public_details <-
  read_delim("extra_data/Public_wgs_info.txt")
annotations_withClinical <-
  annotations_withClinical %>%
  left_join(public_details[, c("Barcode", "Study")], by = c("Barcode")) %>%
  dplyr::mutate(.before = Study.x, Study = if_else(is.na(Study.y), Study.x, Study.y)) %>%
  select(-Study.x, -Study.y)
annotations_withClinical <-
  annotations_withClinical %>%
  left_join(broad_big_list[, c("Barcode", "broad_plate_id")]) %>%
  mutate(batch = if_else(is.na(broad_plate_id), Study, broad_plate_id))

## MEDIAN hSAPIENS COVERAGE DATA
median_coverage <-
  read_delim("extra_data/All_WGS_Median_Coverage.txt") %>%
  select(Barcode, MEDIAN_COVERAGE)
# Grade differentiation annotations
grade_diff_file<-
  read_delim("~/Downloads/FinalHistologySherlock090424.txt")%>%
  filter(!str_detect(Status,"Unresolved"))
grade_diff_info<-
  grade_diff_file%>%
  mutate(Differentiation_result=tolower(Differentiation_result))%>%
  filter(Differentiation_result %in% c("moderate","poorly","well"))%>%
  select(Sherlock_PID,Differentiation_result)
## FULL METADATA
sherlock_metaTable <-
  read_delim("extra_data/sherlock_20240909.csv", col_names = T,delim = ",")%>%
  mutate(STAGE_simple = gsub(pattern = "[ABC].*", replacement = "", x = STAGE_COMPOSITE)) %>%
  mutate(STAGE_simple = if_else(STAGE_simple == 0, "99999", STAGE_simple)) %>%
  mutate(STAGE_simple = if_else(STAGE_simple == "X", "99999", STAGE_simple)) %>%
  mutate(GRADE_DIFFERENTIATION = if_else(GRADE_DIFFERENTIATION == "X", "99999", GRADE_DIFFERENTIATION)) %>%
  mutate(HISTOLOGY_COMPOSITE = if_else(HISTOLOGY_COMPOSITE == "Squamous Cell Carcinoma","Squamous Carcinoma", HISTOLOGY_COMPOSITE)) %>%
  mutate(RACE_DERIVED = case_when(
    !is.na(RACE_DERIVED) ~ RACE_DERIVED,
    STUDY_SITE == "Taiwan" ~ "EAS",
    STUDY_SITE == "Hong Kong" ~ "EAS"
  )) %>%
  mutate(HISTOLOGY_simple = case_when(
    HISTOLOGY_COMPOSITE %in% c("ADENOCARCINOMA", "SQUAMOUS CELL CARCINOMA", "CARCINOID TUMOR") ~ HISTOLOGY_COMPOSITE,
    is.na(HISTOLOGY_COMPOSITE) ~ NA,
    TRUE ~ "Other")) %>%
  `colnames<-`(str_replace(colnames(.), "GENDER", "SEX")) %>%
  `colnames<-`(str_replace(colnames(.), "RACE", "ANCESTRY"))%>%
  mutate(SEX_DERIVED=case_when(SEX_DERIVED==1~"Male",
                               SEX_DERIVED==2~"Female",
                               TRUE~NA)%>%
           as.factor()%>%
           relevel(ref = "Female"))%>%
  mutate(ASTHMA=case_when(ASTHMA==1~"Yes",
                          ASTHMA==2~"No",
                          TRUE~NA)%>%
           as.factor()%>%
           relevel(ref = "No"))%>%
  mutate(PASSIVE_SMOKE=case_when(PASSIVE_SMOKE==1~"Yes",
                                 PASSIVE_SMOKE==2~"No",
                                 TRUE~NA)%>%
           as.factor()%>%
           relevel(ref = "No"))%>%
  mutate(ANY_PREVIOUS_LUNG_DISEASE=case_when(ANY_PREVIOUS_LUNG_DISEASE==1~"Yes",
                                             ANY_PREVIOUS_LUNG_DISEASE==2~"No",TRUE~NA)%>%
           as.factor()%>%
           relevel(ref = "No"))%>%
  left_join(grade_diff_info)%>%
  mutate(GRADE_DIFF_FINAL=relevel(as.factor(Differentiation_result),ref="poorly"))

# excluded samples to remove
sherlock_excluded<-
  read_delim("extra_data/Sherlock_Excluded.csv",delim = ",",col_names = T)%>%
  drop_na(`Sherlock PID`)

wgs_full_annotations <-
  broad_big_list %>%
  filter(Barcode %in% annotations_withClinical$Barcode) %>%
  select(Barcode, sherlock_pid, site_study, tissue_attribute, source_material) %>%
  mutate(
    sherlock_pid = paste0("NSLC-", sherlock_pid),
    Sample_Source = if_else(source_material == "Blood", "Blood", tissue_attribute)
    ) %>%
  left_join(sherlock_metaTable, by = c("sherlock_pid" = "Sherlock_PID"))

pollution_dat <-
  read_delim("extra_data/PollutionData_2023OCT6.csv")
pollution_country <-
  pollution_dat %>%
  summarise(Pollution = median(`Population.Weighted.PM2.5.ug.m3`),.by=c(Study, Country)) %>%
  filter(n() == 1,.by=Study)

s16_metadata <-
  read_delim("~/Downloads/Sherlock_Batch1/NP0493-MB4-manifest.txt") %>%
  separate(`Source PCR Plate`, into = c("PCR_Plate", "Well"), sep = "_(.*?)", extra = "merge") %>%
  mutate_at("PCR_Plate", str_replace_all, "[PC]*", "")
s16_clinical_annos <-
  s16_metadata %>%
  left_join(sherlock_metaTable, by = c("AdditionalAttributes" = "Sherlock_PID")) %>%
  mutate(STAGE_simple = gsub(pattern = "[ABC].*", replacement = "", x = STAGE)) %>%
  mutate(STAGE_simple = if_else(STAGE_simple == 0, "99999", STAGE_simple))

rna_annots <-
  read_delim("~/Sherlock_RNAseq_phase1+2+2R_RNA-seq_clinical_original.txt") %>%
  select(-`Site_REF...18`) %>%
  rename(Site_REF = `Site_REF...3`)

sherlock_rna_batches <-
  read_delim("extra_data/simple_sherlock_rnaBatch.txt", col_names = F) %>%
  rename(Date = "X1", Machine = "X2", SampleID = "X3") %>%
  mutate(batch = paste0(Date, "_", Machine)) %>%
  group_by(SampleID) %>%
  filter(n() < 2 | Date == 210125) %>%
  select(RNAseq_SampleID = "SampleID", batch) %>%
  mutate_at("batch", as.character) %>%
  ungroup()
tcga_rna_batches <-
  rna_annots %>%
  select(RNAseq_SampleID) %>%
  filter(str_detect(RNAseq_SampleID, "TCGA")) %>%
  mutate(col2 = str_split(RNAseq_SampleID, pattern = "-")) %>%
  unnest_wider(col2, names_sep = "minimal") %>%
  select(1, 3, 7, 8) %>%
  `colnames<-`(c("RNAseq_SampleID", "Hosp", "Plate", "Center")) %>%
  mutate(batch = paste0(Center, "_", Plate)) %>%
  select(RNAseq_SampleID, batch)
tcga_rna_batchInfo <-
  rna_annots %>%
  select(RNAseq_SampleID) %>%
  filter(str_detect(RNAseq_SampleID, "TCGA")) %>%
  mutate(col2 = str_split(RNAseq_SampleID, pattern = "-")) %>%
  unnest_wider(col2, names_sep = "minimal") %>%
  select(1, 3, 7, 8) %>%
  `colnames<-`(c("RNAseq_SampleID", "Hosp", "Plate", "Center"))
EAGLE_batch_info <-
  rna_annots %>%
  filter(Dataset == "EAGLE") %>%
  select(RNAseq_SampleID, FastqFile) %>%
  mutate(batch = str_split(FastqFile, "/") %>%
           lapply(FUN = function(x) {
             x[length(x) - 1]
           }) %>% unlist()) %>%
  select(RNAseq_SampleID, batch)
rna_batch_info <-
  bind_rows(sherlock_rna_batches, EAGLE_batch_info, tcga_rna_batches)
rna_annots <-
  left_join(rna_annots, rna_batch_info) %>%
  mutate(batch = if_else(is.na(batch), Dataset, batch))
rna_annots_full <-
  rna_annots %>%
  select(RNAseq_SampleID:Dataset_NUM, batch) %>%
  left_join(sherlock_metaTable) %>%
  mutate(Site_REF = str_replace(Site_REF, pattern = "CCR-Anish", replacement = "NCI-CCR"))
rna_annots_full%>%
  left_join(read_delim("extra_data/RNA_extraction_batch.txt")%>%
            select(`Sherlock PID`,`Tissue Attribute`,`Created By`,`Source Material`)%>%
            unique(),
          by=c("Sherlock_PID"="Sherlock PID","Type"="Tissue Attribute"))

conversion_guide <-
  left_join(
    s16_metadata %>%
      select(SampleID, AdditionalAttributes, `Tumor-NormalStatus`) %>%
      mutate(AdditionalAttributes = gsub("NSLC-", "", AdditionalAttributes)),
    broad_big_list,
    by = c("AdditionalAttributes" = "sherlock_pid", "Tumor-NormalStatus" = "tissue_attribute")
  ) %>%
  select(SampleID, source_material, Barcode) %>%
  filter(source_material != "Blood") %>%
  select(-source_material)

################################################################################
# READ KRAKEN OUTPUTS
################################################################################
kraken_raw <-
  read_delim("/Volumes/Sherlock_Lung/JohnMce/Kraken_WGS3/nov24_test_collated.txt") %>%
  rename(tax_id = `4`) %>%
  select(tax_id, any_of(annotations_withClinical$Barcode)) %>%
  mutate_at("tax_id", as.character)
# Remove public data; too noisy
wgs_nopublic<-
  annotations_withClinical%>%
  filter(Smoking=="Non-Smoker",
         !Study %in%c("TCGA-legacy","TCGA-new","EAGLE_Smoker","EAGLE_Nonsmoker",
                      "EAGLE_LUSC_LUNE","Lee-2017","Lee-2019","Leong-2019",
                      "Imielinski"))%>%
  left_join(wgs_full_annotations,by=c("Barcode"))%>%
  filter(!str_detect(HISTOLOGY_COMPOSITE,"EXCLUDE"))%>%
  filter(!sherlock_pid %in% sherlock_excluded$`Sherlock PID`)
wgs_sherlockOnly<-
  kraken_raw%>%
  select(any_of(c("tax_id",wgs_nopublic$Barcode)))

# Read RNA Data
rna <-
  read_delim("/Volumes/Sherlock_Lung/JohnMce/Kraken_RNA3/rna_nov24_test_collated.txt")%>%
  rename(tax_id = "4")%>%
  mutate_at("tax_id", as.character)

# filter to only fresh-frozen samples with RNA extracted at nationwide, remove all identified metastases
rna_nationwide_frozen_Samples <-
  rna_annots%>%
  left_join(read_delim("extra_data/RNA_extraction_batch.txt")%>%
              select(`Sherlock PID`,`Tissue Attribute`,`Created By`,`Source Material`)%>%
              unique(),
            by=c("Sherlock_PID"="Sherlock PID","Type"="Tissue Attribute"))%>%
  filter(str_detect(Dataset,"Sherlock"))%>%
  drop_na(`Sherlock_PID`)%>%
  replace_na(list(Comments="none"))%>%
  filter(`Created By`=="Nationwide",
         `Source Material`=="Frozen Tissue",
         Comments %in% c("none","Tumor Only"),
         str_detect(HISTOLOGY_COMPOSITE,"EXCLUDE", negate=TRUE),
         !Sherlock_PID %in% sherlock_excluded$`Sherlock PID`)
rna_NW_FF<-
  rna%>%
  select(any_of(c("tax_id",rna_nationwide_frozen_Samples$RNAseq_SampleID)))

# 16s kraken
s16_kraken <-
  #read_delim("/Volumes/Sherlock_Lung/JohnMce/Kraken_16S2/s16_ncbi_bracken0.02_collated.txt")
  read_delim("/Volumes/Sherlock_Lung/JohnMce/Kraken_16S3/s16_nov24_test_collated.txt") %>%
  rename(tax_id = "4") %>%
  replace(is.na(.), 0) %>%
  #filter(tax_id %in% kraken_microbeTaxonomy$tax_id)%>%
  `colnames<-`(str_replace(colnames(.), "Sample_", "") %>% str_replace("-[TCGA]*$", "")) %>%
  mutate_at("tax_id", as.character)
s16_to_wgs <-
  read_delim("~/Documents/Sherlock_16s_withPIDs.txt") %>%
  filter(`Sample ID` %in% s16_metadata$SampleID) %>%
  select(`Sherlock PID`, `Sample ID`, `Tumor\n`) %>%
  mutate(`Sherlock PID` = str_replace(`Sherlock PID`, "NSLC-", "")) %>%
  left_join(
    broad_big_list[, c("sherlock_pid", "tissue_attribute", "source_material", "Barcode")] %>%
      mutate(tissue_attribute = gsub("Peritumoral", "Tumor", tissue_attribute)),
    by = c("Sherlock PID" = "sherlock_pid", "Tumor\n" = "tissue_attribute")
  ) %>%
  filter(source_material == "Lung")

## Not all mutational signatures were included in this table so I added them back
## maybe some are unreliable?
sherlock_data_full <-
  sherlock_sig_cosmic %>%
  pivot_longer(-1, names_to = "Gene", values_to = "Alteration") %>%
  mutate(Alteration = if_else(Alteration > 0, "Yes", "No"), Type = "Signature_Cosmic") %>%
  left_join(annotations_withClinical[, c("Barcode", "Subject")], by = c("Tumor_Barcode" = "Barcode")) %>%
  select(colnames(sherlock_data_full)) %>%
  filter(!Gene %in% sherlock_data_full$Gene) %>%
  rbind(sherlock_data_full)

# Sample read depth info
public_read_depth<-read_delim("extra_data/public_dat_readsums.txt")
read_depth_wgs<-read_delim("extra_data/AlignmentStats.txt")
wgs_read_depth<-
  read_depth_wgs%>%
  select(Barcode:TOTAL_READS)%>%
  rbind(public_read_depth%>%
          rename(Barcode="Sample")%>%
          select(Barcode:TOTAL_READS))%>%filter(CATEGORY=="FIRST_OF_PAIR")

### Database information
salter_list <-
  read_delim("extra_data/Salter_blacklist.txt", delim = "\t") %>%
  left_join(kraken_microbeTaxonomy, by = c("Genus" = "name")) %>%
  filter(type == "genus")
# GRIMER_db <-
#   read_delim("extra_data/GRIMER_db_ncbi.txt",delim = ", ",col_names = F)%>%
#   `colnames<-`("tax_id")
pathogen_v_host <-
  read_delim("extra_data/PathogenVsHostDB-2019-05-30.csv")
salter_list_nonHuman <-
  salter_list %>%
  mutate(Human_associated = (Genus %in% (pathogen_v_host %>%
                                           filter(HostSpecies == "Homo sapiens") %>%
                                           count(Genus) %>%
                                           filter(n > 1) %>%
                                           pull(Genus)))) %>%
  filter(!Human_associated)