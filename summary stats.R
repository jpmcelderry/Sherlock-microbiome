data_upset_plot%>%
  filter(Type=="Tumor")%>%
  count(SEX_DERIVED)
data_upset_plot%>%
  filter(Type=="Tumor")%>%
  summarise(median(AGE_AT_DIAGNOSIS,na.rm=TRUE))
data_upset_plot%>%
  filter(Type=="Tumor")%>%
  count(str_detect(Data,"WGS"),ANCESTRY_DERIVED)
data_upset_plot%>%
  filter(Type=="Tumor")%>%
  count(STUDY_SITE)%>%
  arrange(-n)
data_upset_plot%>%
  filter(Type=="Tumor")%>%
  mutate(HISTOLOGY_COMPOSITE=fct_lump(HISTOLOGY_COMPOSITE,n = 3))%>%
  count(HISTOLOGY_COMPOSITE)%>%
  arrange(-n)
data_upset_plot%>%
  mutate(genomics=paste0(Data,collapse = ";"))%>%
  count(Type)
data_upset_plot%>%
  mutate(genomics=paste0(Data,collapse = ";"))%>%
  count(Type,Data)

kraken_assignment_results<-
  RNA_kingdoms%>%
  full_join(WGS_kingdoms)%>%
  full_join(s16_kingdoms)%>%
  replace(is.na(.),0)%>%
  pivot_longer(-1,names_to = "Barcode",values_to = "Reads")%>%
  mutate(platform=case_when(Barcode %in% colnames(s16_scrubbed)~"16S",
                            Barcode %in% colnames(rna_combatd_decontamd)~"RNA",
                            Barcode %in% colnames(wgs_raw_combined)~"WGS"))%>%
  left_join(data_upset_plot%>%
              unnest_longer(IDs)%>%
              select(Barcode="IDs",Type))%>%
  filter(Type %in% c("Tumor","Normal"),!tax_id=="9606")%>%
  mutate(Reads=Reads/sum(Reads),.by = c(Barcode))%>%
  summarise(median_pct=median(Reads)*100,.by = c(platform,tax_id))%>%
  left_join(unified_taxonomy)
kraken_assignment_results

kraken_assignment_results_microbial<-
  RNA_kingdoms%>%
  full_join(WGS_kingdoms)%>%
  full_join(s16_kingdoms)%>%
  replace(is.na(.),0)%>%
  pivot_longer(-1,names_to = "Barcode",values_to = "Reads")%>%
  mutate(platform=case_when(Barcode %in% colnames(s16_scrubbed)~"16S",
                            Barcode %in% colnames(rna_combatd_decontamd)~"RNA",
                            Barcode %in% colnames(wgs_raw_combined)~"WGS"))%>%
  left_join(data_upset_plot%>%
              unnest_longer(IDs)%>%
              select(Barcode="IDs",Type))%>%
  filter(Type %in% c("Tumor","Normal"),!tax_id %in% c("0","2759","9605"))%>%
  mutate(Reads=Reads/sum(Reads),.by = c(Barcode))%>%
  summarise(median_pct=median(Reads)*100,.by = c(platform,tax_id))%>%
  left_join(unified_taxonomy)
kraken_assignment_results_microbial

total_reads<-
  rna_rawCounts_decontamd%>%
  full_join(s16_scrubbed)%>%
  full_join(wgs_raw_combined)%>%
  left_join(unified_taxonomy)%>%
  filter(type=="superkingdom",str_detect(taxonomy,"Bacteria"))%>%
  select(where(is.numeric))%>%
  replace(is.na(.),0)%>%
  colSums()%>%
  as.data.frame()%>%
  rownames_to_column("Barcode")%>%rename(Reads=".")%>%
  filter(Barcode %in% c(colnames(s16_scrubbed),
                        colnames(rna_combatd_decontamd),
                        colnames(wgs_raw_combined)))%>%
  left_join(all_sample_TNstatus)%>%
  group_by(dataset,Source)%>%
  summarise(median(Reads),n=n())
total_reads

relative_abundance_summary<-
  s16_scrubbed%>%
  left_join(rna_rawCounts_decontamd)%>%
  left_join(wgs_raw_combined)%>%
  replace(is.na(.),0)%>%
  pivot_longer(where(is.numeric),names_to = "Barcode",values_to = "Reads")%>%
  left_join(unified_taxonomy)%>%
  mutate(Reads=Reads/sum(Reads),.by = c(Barcode,type))%>%
  drop_na(Reads)%>%
  filter(type %in% c("superkingdom","phylum","genus"))%>%
  mutate(platform=case_when(Barcode %in% colnames(s16_scrubbed)~"16S",
                            Barcode %in% colnames(rna_combatd_decontamd)~"RNA",
                            Barcode %in% colnames(wgs_raw_combined)~"WGS"))
relative_abundance_summary

phylum_relative_abundance_ranges<-
  relative_abundance_summary%>%
  filter(name %in% c("Pseudomonadota","Actinomycetota","Bacillota"))%>%
  left_join(data_upset_plot%>%
              unnest_longer(IDs)%>%
              select(Barcode="IDs",Type))%>%
  summarise(mean=mean(Reads)*100,
            median=median(Reads)*100,
            min=min(Reads)*100,
            max=max(Reads)*100,
            .by = c(platform,name))
phylum_relative_abundance_ranges

big_bacteria_relative_abundances<-
  relative_abundance_summary%>%
  filter(name %in% c("Acinetobacter","Corynebacterium","Pseudomonas","Staphylococcus","Streptococcus"))%>%
  left_join(data_upset_plot%>%
              unnest_longer(IDs)%>%
              select(Barcode="IDs",Type))%>%
  filter(Type %in% c("Tumor","Normal"))%>%
  summarise(mean=mean(Reads)*100,
            median=median(Reads)*100,
            min=min(Reads)*100,
            max=max(Reads)*100,
            .by = c(platform,name))
big_bacteria_relative_abundances
big_bacteria_relative_abundances%>%
  summarise(round(range(mean),digits = 1),.by=name)

BAL_bacteria_relative_abundances<-
  relative_abundance_summary%>%
  filter(name %in% c("Prevotella","Haemophilus","Veillonella","Neisseria","Streptococcus"))%>%
  left_join(data_upset_plot%>%
              unnest_longer(IDs)%>%
              select(Barcode="IDs",Type))%>%
  group_by(platform,Barcode)%>%
  summarise(sum=sum(Reads))%>%
  summarise(mean=mean(sum)*100)
BAL_bacteria_relative_abundances

rm(relative_abundance_summary)