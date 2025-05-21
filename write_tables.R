# Write sample metadata
cohort_meta_table<-
  data_upset_plot%>%
  mutate(Data=sapply(Data,function(x)paste0(x,collapse = ", ")))%>%
  mutate(IDs=sapply(IDs,function(x)paste0(x,collapse = ", ")))%>%
  select(Sherlock_PID,
         Type,Data,IDs,
         STUDY_SITE,
         HISTOLOGY_COMPOSITE,
         STAGE_simple,
         AGE_AT_DIAGNOSIS,
         SEX_DERIVED,
         ANCESTRY_DERIVED,
         PASSIVE_SMOKE,
         SURVIVAL_STATUS="VITAL_STATUS",
         SURVIVAL_TIME_WEEKS_DERIVED,
         ANY_PREVIOUS_LUNG_DISEASE,
         RECURRENCE,
         METASTASIS)%>%
  mutate(SURVIVAL_STATUS=case_when(SURVIVAL_STATUS==1~"Deceased",
                                   SURVIVAL_STATUS==0~"Alive",
                                   TRUE ~ NA))%>%
  mutate(RECURRENCE=case_when(RECURRENCE==1~"Yes",
                              RECURRENCE==2~"No",
                              TRUE ~ NA))%>%
  mutate(STAGE_simple=if_else(STAGE_simple=="99999",NA,STAGE_simple))%>%
  mutate(STUDY_SITE=case_when(STUDY_SITE=="Mayo"~"Minnesota, USA",
                           STUDY_SITE=="Roswell Park"~"New York, USA",
                           STUDY_SITE=="INCan"~"Mexico City, Mexico",
                           STUDY_SITE=="Yale"~"Connecticut, USA",
                           STUDY_SITE=="Harvard"~"Massachusetts, USA",
                           STUDY_SITE=="Moffitt"~"Florida, USA",
                           STUDY_SITE=="Laval"~"Quebec, Canada",
                           STUDY_SITE=="IARC"~"IARC (Serbia, Czech Republic, Romania, Poland, Russia)",
                           STUDY_SITE=="Valencia"~"Valencia, Spain",
                           STUDY_SITE=="Peru"~"Lima, Peru",
                           STUDY_SITE=="Toronto"~"Toronto, Canada",
                           STUDY_SITE=="Nice"~"Nice, France",
                          .default=STUDY_SITE)
         )%>%
  mutate(ANCESTRY_DERIVED=case_when(ANCESTRY_DERIVED=="EUR"~"European (EUR)",
                                    ANCESTRY_DERIVED=="AFR"~"African (AFR)",
                                    ANCESTRY_DERIVED=="EAS"~"East Asian (EAS)",
                                    ANCESTRY_DERIVED=="AMR or Mixed"~"American (AMR) or Mixed",
                                    TRUE ~ ANCESTRY_DERIVED)
         )%>%
  mutate_at("HISTOLOGY_COMPOSITE",str_to_sentence)

colnames(cohort_meta_table)<-
  str_replace_all(colnames(cohort_meta_table),"_"," ")%>%
  str_to_title()%>%
  str_replace(" Derived","")%>%
  str_replace(" Simple","")

cohort_table_1_columns<-
  c("Sex",
    "Ancestry",
    "Study Site",
    "Stage",
    "Histology Composite")

# Main table 1
cohort_meta_table%>%
  separate_longer_delim(Data,", ")%>%
  mutate(`Paired normal`="Normal" %in% Type,
         `Paired blood`="Blood" %in% Type,
         Unpaired = !(`Paired normal`|`Paired blood`),
         .by=c(`Sherlock Pid`,Data))%>%
  filter(Type=="Tumor")%>%
  mutate(`Histology Composite`=fct_lump_n(`Histology Composite`,n = 4))%>%
  tbl_summary(by = Data,
              type = list(cohort_table_1_columns ~ "categorical"),
              include = c(`Age At Diagnosis`,cohort_table_1_columns)
  )%>%
  bold_labels()%>%
  as_gt()%>%
  gt::gtsave("Figure_1_table.docx","plots/current_plots/")

write_tsv(x = cohort_meta_table,"supplement_tables/metadata_suppleTable.tsv")

# Write raw counts tables
s16_kraken%>%
  select_if(str_detect(colnames(.),"NTC")|str_detect(colnames(.),"Water")|colnames(.) %in% colnames(s16_scrubbed))%>%
  replace(.==0,NA)%>%
  right_join(unified_taxonomy,.)%>%
  write_tsv("supplement_tables/s16_raw.tsv",na = "")

kraken_raw%>%
  select(colnames(wgs_decontamd))%>%
  replace(.==0,NA)%>%
  right_join(unified_taxonomy,.)%>%
  write_tsv("supplement_tables/wgs_raw.tsv",na = "")
  
rna%>%
  select(colnames(rna_combatd_decontamd))%>%
  replace(.==0,NA)%>%
  right_join(unified_taxonomy,.)%>%
  filter(!str_detect(type,"species"))%>%
  write_tsv("supplement_tables/rna_raw.tsv",na = "")

## cleaned data
s16_scrubbed%>%
  right_join(unified_taxonomy,.)%>%
  replace(.==0,NA)%>%
  write_tsv("supplement_tables/s16_cleaned.tsv",na = "")
wgs_decontamd%>%
  rbind(wgs_decontamd_phylum)%>%
  full_join(rbind(NG232_decontamd,
                  NG232_decontamd_phylum))%>%
  right_join(unified_taxonomy,.)%>%
  replace(.==0,NA)%>%
  write_tsv("supplement_tables/wgs_cleaned.tsv",na="")
rna_combatd_decontamd%>%
  rbind(rna_combatd_decontamd_phylum)%>%
  replace(.==0,NA)%>%
  right_join(unified_taxonomy,.)%>%
  filter(!str_detect(type,"species"))%>%
  write_tsv("supplement_tables/rna_cleaned.tsv",na="")

# ICC results
ICC_beta_WGS_16s%>%tidy()%>%
  filter(term=="Sherlock_PID")%>%
  mutate()%>%
  mutate(Experiment="WGS v 16S")%>%
  rbind(ICC_beta_RNA_16s%>%
          tidy()%>%
          filter(term=="Sherlock_PID")%>%
          mutate(Experiment="RNA v 16S"))%>%
  rbind(ICC_beta_RNA_WGS%>%
          tidy()%>%
          filter(term=="Sherlock_PID")%>%
          mutate(Experiment="RNA v WGS"))%>%
  rbind(ICC_beta_AllThree%>%
          tidy()%>%
          filter(term=="Patient_tissue")%>%
          mutate(Experiment="RNA v WGS v 16S"))%>%
  relocate(Experiment,.before = term)%>%
  write_tsv("supplement_tables/betadiversity_ICC.tsv")

## Contaminant list
salter_list_nonHuman%>%
  write_tsv("supplement_tables/contaminators_list.tsv")

