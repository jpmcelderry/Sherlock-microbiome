# Write sample metadata
data_upset_plot%>%
  mutate(Data=sapply(Data,function(x)paste0(x,collapse = ", ")))%>%
  mutate(IDs=sapply(IDs,function(x)paste0(x,collapse = ", ")))%>%
  select(Sherlock_PID,Type,Data,IDs,STUDY_SITE,HISTOLOGY_COMPOSITE,
         STAGE_DERIVED,GRADE_DIFF_FINAL,AGE_AT_DIAGNOSIS,
         SEX_DERIVED,ANCESTRY_DERIVED,PASSIVE_SMOKE,VITAL_STATUS,
         SURVIVAL_TIME_WEEKS_DERIVED,ANY_PREVIOUS_LUNG_DISEASE,RECURRENCE,ASTHMA)%>%
  write_tsv("supplement_tables/metadata_suppleTable.tsv")

# Write raw counts tables
s16_kraken%>%
  select_if(str_detect(colnames(.),"NTC")|str_detect(colnames(.),"Water")|colnames(.) %in% colnames(s16_scrubbed))%>%
  replace(.==0,NA)%>%
  right_join(taxonomy_16S,.)%>%
  write_tsv("supplement_tables/s16_raw.tsv")

kraken_raw%>%
  select(colnames(wgs_decontamd))%>%
  replace(.==0,NA)%>%
  right_join(taxonomy_WGS,.)%>%
  write_tsv("supplement_tables/wgs_raw.tsv")
  
rna%>%
  select(colnames(rna_combatd_decontamd))%>%
  replace(.==0,NA)%>%
  right_join(taxonomy_RNA,.)%>%
  filter(!str_detect(type,"species"))%>%
  write_tsv("supplement_tables/rna_raw.tsv")

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
salter_list%>%
  select(Genus)%>%
  left_join(kraken_taxonomy,by=c("Genus"="name"))%>%
  rename(Genus="name")%>%
  left_join(pathogen_v_host%>%
              filter(HostSpecies=="Homo sapiens")%>%
              count(Genus))%>%
  rename(n_species="n")%>%
  replace(is.na(.),0)%>%
  mutate(Human_associated=n_species>1)%>%
  write_tsv("supplement_tables/contaminators_list.tsv")

