######## FIGURE 1
## 1a
wgs_raw_combined<-
  wgs_rawCounts_decontamd%>%
  full_join(ng232_rawCounts_decontamd)%>%
  replace(is.na(.),0)

data_upset_plot <-
  rna_annots_full %>%
  filter(RNAseq_SampleID %in% c(colnames(rna_combatd_decontamd))) %>%
  select(Sherlock_PID, Type, ANCESTRY_DERIVED, STUDY_SITE,Barcode="RNAseq_SampleID") %>%
  mutate(Data = "RNA") %>%
  distinct(Sherlock_PID,Type,.keep_all = T) %>%
  rbind(s16_clinical_annos %>%
          filter(SampleID %in% colnames(s16_scrubbed)) %>%
          select(Sherlock_PID = "AdditionalAttributes", Type = "Tumor-NormalStatus", ANCESTRY_DERIVED, STUDY_SITE,Barcode="SampleID") %>%
          mutate(Data = "16S") %>%
          distinct(Sherlock_PID,Type,.keep_all = T)) %>%
  rbind(wgs_full_annotations %>%
          filter(Barcode %in% colnames(wgs_raw_combined)) %>%
          select(Sherlock_PID = "sherlock_pid", Type = "Sample_Source", ANCESTRY_DERIVED, STUDY_SITE, Barcode) %>%
          mutate(Data = "WGS") %>%
          distinct(Sherlock_PID,Type,.keep_all = T)) %>%
  mutate(Type = str_replace(Type, "Lung", "Normal")) %>%
  summarise(Data = list(Data),IDs=list(Barcode),.by=c(Sherlock_PID, Type)) %>%
  left_join(sherlock_metaTable) %>%
  filter(Type %in% c("Tumor", "Normal", "Blood")) %>%
  rowwise() %>%
  group_by(Type, Data) %>%
  mutate(n = n()) %>%
  ungroup()

sankey_in <-
  data_upset_plot %>%
  drop_na(Type, Data) %>%
  mutate(sample_status = sapply(Data, 
                                function(x){
                                  paste0(x, collapse = " +\n")
                                  }
                                )) %>%
  mutate(Overlap = if_else(str_detect(sample_status, "\\+"),
                                   sample_status,
                                   paste0(sample_status, " only")
  )) %>%
  mutate_at("Overlap", fct_lump_n, n = 4, other_level = "All other") %>%
  unnest_longer(col = Data)%>%
  make_long(Data, `Tissue Type` = "Type",Overlap) %>%
  mutate_at("next_node", factor, levels = c("All other",
                                       "WGS only",
                                       "16S +\nWGS",
                                       "RNA +\n16S",
                                       "RNA +\n16S +\nWGS",
                                       "Blood",
                                       "Normal",
                                       "Tumor",
                                       "WGS",
                                       "RNA",
                                       "16S"))%>%
  mutate_at("node", factor, levels = c("All other",
                                       "WGS only",
                                       "16S +\nWGS",
                                       "RNA +\n16S",
                                       "RNA +\n16S +\nWGS",
                                       "Blood",
                                       "Normal",
                                       "Tumor",
                                       "WGS",
                                       "RNA",
                                       "16S"))

sankey_plot <-
  ggplot(sankey_in, aes(
    x = x,
    next_x = next_x,
    node = node,
    next_node = next_node,
    fill = factor(node)
  )) +
  geom_sankey(
    flow.alpha = .6,
    node.color = "gray30"
  ) +
  geom_sankey_label(
    aes(x = as.numeric(x),
        label = after_stat(paste0(node, "\n(n = ", freq, ")"))),
    size = 2, color = "white", fill = "gray40"
  ) +
  scale_fill_viridis_d(drop = FALSE) +
  theme_sankey(base_size = 12, base_family = "Roboto Condensed") +
  theme(legend.position = "none") +
  labs(x = NULL) +
  scale_x_discrete(expand = c(0.2, 0),position = "top") +
  theme(legend.position = "none")
sankey_plot

### Fig 1b 
workflow_grob <-
  grid::rasterGrob(png::readPNG(source = "plots/Microbiome_workflow_feb25.png"),
                   interpolate = FALSE, gp = gpar(lwd = 2, col = "white", fill = "#00000000")
  )
workflow <- ggplot() +
  geom_blank() +
  annotation_custom(workflow_grob,
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  theme(plot.background = element_rect(color = "white")) +
  theme_classic()
workflow

## Fig 1c
highest_level_clades<-
  kraken_taxonomy%>%
  filter(type %in% c("kingdom","superkingdom")|name=="Homo")%>%
  pull(tax_id)
WGS_kingdoms<-
  read_delim("/Volumes/Sherlock_Lung/JohnMce/Kraken_WGS3/kraken0.1-collated.txt")%>%
  rename(tax_id=`4`)%>%
  filter(tax_id %in% highest_level_clades|tax_id==0)%>%
  select(any_of(colnames(wgs_raw_combined)))%>%
  select(-any_of(annotations_withClinical%>%
                   filter(!Sample_Source %in% c("Lung","Tumor"))%>%
                   pull(Barcode)))%>%
  mutate_at("tax_id",as.character)
RNA_kingdoms<-
  rna%>%
  filter(tax_id %in% highest_level_clades|tax_id==0)%>%
  select(any_of(colnames(rna_combatd_decontamd)))%>%
  mutate_at("tax_id",as.character)
s16_kingdoms<-
  read_delim("/Volumes/Sherlock_Lung/JohnMce/Kraken_16S3/kraken0.02-collated.txt")%>%
  rename(tax_id=`4`)%>%
  `colnames<-`(str_replace(colnames(.),"Sample_","")%>%
                 str_replace("-[TCGA]*$",""))%>%
  filter(tax_id %in% highest_level_clades|tax_id==0)%>%
  select(any_of(colnames(s16_scrubbed)))%>%
  mutate_at("tax_id",as.character)
#### Bar plots by kingdom
kingdoms_bar<-
  RNA_kingdoms%>%
  full_join(WGS_kingdoms)%>%
  full_join(s16_kingdoms)%>%
  replace(is.na(.),0)%>%
  filter(rowMeans(.>0,)>0.1)%>%
  pivot_longer(2:ncol(.),names_to = "Barcode",values_to = "Reads")%>%
  mutate(Experiment=case_when(Barcode %in% colnames(RNA_kingdoms)~"RNA-seq",
                              Barcode %in% colnames(WGS_kingdoms)~"WGS",
                              Barcode %in% colnames(s16_kingdoms)~"16S"))%>%
  left_join(kraken_taxonomy)

kingdom_barplot<-
  kingdoms_bar%>%
  replace_na(list(name="unknown"))%>%
  filter(!name %in% c("Metazoa","Sangervirae","Orthornavirae","Eukaryota"))%>%
  ggplot(aes(x=fct_relevel(name,"unknown",after = Inf),y=log10(Reads+1),fill=Experiment))+
  geom_boxplot(outlier.alpha = 0.2)+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1))+
  scale_fill_brewer(type = "qual")+
  labs(x="",fill="dataset",y="log10 Median Reads")
kingdom_barplot
# ggsave("plots/current_plots/Read_totals.pdf",device=cairo_pdf,width=8,height=4)

###### BACTERIAL PROPORTIONS - Fig 1d
wgs_bacterial_proportions<-
  read_delim("extra_data/AlignmentStats.txt")%>%
  filter(CATEGORY=="PAIR")%>%
  select(Barcode,TOTAL_READS)%>%
  left_join(annotations_withClinical)%>%
  filter(Barcode %in% colnames(wgs_raw_combined))%>%
  left_join(kraken_raw%>%
              filter(tax_id=="2")%>%
              pivot_longer(-1,names_to = "Barcode",values_to = "Bacteria"))%>%
  drop_na(Bacteria)%>%
  mutate(Bact_proportion=Bacteria/TOTAL_READS*1e6)

rna_total_reads<-
  read_delim("extra_data/RNAseq_QC.txt")%>%
  mutate_at("Features",trimws)%>%
  filter(Features=="Number of input reads")%>%
  pivot_longer(-Features,names_to = "RNAseq_SampleID",values_to = "TOTAL_READS")%>%
  select(-Features)

rna_bacterial_proportions<-
  rna%>%
  filter(tax_id=="2")%>%
  pivot_longer(-1,names_to = "RNAseq_SampleID",values_to = "Bacteria")%>%
  left_join(rna_total_reads)%>%
  mutate(Bact_proportion=Bacteria/as.numeric(TOTAL_READS)*1e6)%>%
  left_join(rna_annots)

s16_total_reads<-
  read_delim("extra_data/read_counts_16S.txt2")%>%
  mutate_at("Sample",str_replace,"Sample_","")%>%
  mutate_at("Sample",str_replace,"-[ATCG]*$","")%>%
  rename(SampleID="Sample")

s16_bacterial_proportions<-
  s16_total_reads%>%
  left_join(s16_kraken%>%
              filter(tax_id=="2")%>%
              column_to_rownames("tax_id")%>%
              t()%>%
              as.data.frame()%>%
              rownames_to_column("SampleID"))%>%
  left_join(s16_metadata)%>%
  drop_na(`2`)%>%
  mutate(Bact_proportion=`2`/TOTAL_READS*1e6)

bacterial_proportions_data<-
  wgs_bacterial_proportions%>%
  select(Barcode,Sample_Source,Bact_proportion)%>%
  mutate(Experiment="WGS")%>%
  mutate_at("Sample_Source",str_replace,"Lung","Normal")%>%
  rbind(rna_bacterial_proportions%>%
          select(Barcode="RNAseq_SampleID",Sample_Source="Type",Bact_proportion)%>%
          mutate(Experiment="RNA-seq")%>%
          filter(Barcode %in% colnames(rna_combatd_decontamd)))%>%
  rbind(s16_bacterial_proportions%>%
          select(Barcode="SampleID",Sample_Source="Tumor-NormalStatus",Bact_proportion)%>%
          mutate(Experiment="16S")%>%
          filter(Barcode %in% colnames(s16_scrubbed)))%>%
  filter(Sample_Source %in% c("Blood","Normal","Tumor"))

bacterial_proportions_boxplot<-
  bacterial_proportions_data%>%
  ggplot(aes(x=factor(Sample_Source,levels=c("Tumor","Normal","Blood")),
             y=log10(Bact_proportion),fill=Experiment))+
  geom_boxplot(position=position_dodge2(preserve = "single"),width=0.6,outliers = F)+
  labs(x=NULL,y="log10 Bacterial Reads Per Million")+
  scale_fill_brewer(type = "qual")+
  theme(axis.title.y = element_text(size=10))
bacterial_proportions_boxplot

## 1e
## Compile samples that are shared between platforms
paired_rna_wgs<-
  rna_annots%>%
  mutate(Sherlock_PID=paste0(Sherlock_PID,"_",Type))%>%
  select(SampleID=RNAseq_SampleID,Sherlock_PID)%>%
  mutate(Sherlock_PID=str_replace(Sherlock_PID,"NSLC-",""),experiment="RNA")%>%
  rbind(broad_big_list%>%
          filter(source_material=="Lung")%>%
          mutate(Sherlock_PID=paste0(sherlock_pid,"_",tissue_attribute))%>%
          select(SampleID=Barcode,Sherlock_PID)%>%
          mutate(experiment="WGS"))

paired_wgs_16s<-
  broad_big_list%>%
  filter(source_material=="Lung")%>%
  mutate(Sherlock_PID=paste0(sherlock_pid,"_",tissue_attribute))%>%
  select(SampleID=Barcode,Sherlock_PID)%>%
  mutate(experiment="WGS")%>%
  rbind(s16_metadata%>%
          mutate(Sherlock_PID=paste0(AdditionalAttributes,"_",`Tumor-NormalStatus`))%>%
          select(SampleID,Sherlock_PID)%>%
          mutate(Sherlock_PID=str_replace(Sherlock_PID,"NSLC-",""),experiment="16s"))

samples_crossPlatform<-
  paired_wgs_16s%>%
  pivot_wider(names_from = "experiment",values_from = "SampleID")%>%
  full_join(paired_rna_wgs%>%
              pivot_wider(names_from = "experiment",values_from = "SampleID"))%>%
  drop_na()%>%
  mutate_at(c("WGS","16s","RNA"),as.character)%>%
  filter(WGS %in% colnames(wgs_raw_combined),
         `16s` %in% colnames(s16_scrubbed),
         RNA %in% colnames(rna_combatd_decontamd))%>%
  pivot_longer(c("WGS","16s","RNA"),names_to = "platform",values_to = "Barcode")

### Create table with all tumor normal status info, cross-platform
all_sample_TNstatus<-
  c(annotations_withClinical$Barcode,
    s16_metadata$SampleID,
    rna_annots$RNAseq_SampleID)%>%
  as.data.frame()%>%
  rename(Barcode=".")%>%
  left_join(annotations_withClinical[,c("Barcode","Sample_Source")]%>%
              mutate(type="WGS"))%>%
  mutate(Sample_Source=if_else(Sample_Source=="Lung","Normal",Sample_Source))%>%
  left_join(rna_annots[,c("RNAseq_SampleID","Type")]%>%
              mutate(type="RNA"),by=c("Barcode"="RNAseq_SampleID"))%>%
  mutate(Source=if_else(is.na(Sample_Source),Type,Sample_Source),
         type=if_else(is.na(type.x),type.y,type.x))%>%
  select(Barcode,type,Source)%>%
  left_join(s16_metadata[,c("SampleID","Tumor-NormalStatus")]%>%
              mutate(type="16s")%>%
              replace_na(list(`Tumor-NormalStatus`="NTC")),
            by=c("Barcode"="SampleID"))%>%
  mutate(Source=if_else(is.na(Source),`Tumor-NormalStatus`,Source),
         type=if_else(is.na(type.x),type.y,type.x))%>%
  select(Barcode,dataset=type,Source)

compare_read_counts<-
  rna%>%
  full_join(s16_kraken)%>%
  full_join(kraken_raw)%>%
  left_join(kraken_taxonomy)%>%
  filter(type=="superkingdom",str_detect(taxonomy,"Bacteria"))%>%
  select(where(is.numeric))%>%
  replace(is.na(.),0)%>%
  colSums()%>%
  as.data.frame()%>%
  rownames_to_column("Barcode")%>%rename(Reads=".")%>%
  filter(Barcode %in% c(colnames(s16_scrubbed),
                        colnames(rna_combatd_decontamd),
                        colnames(wgs_raw_combined)))%>%
  left_join(all_sample_TNstatus)
# plot read totals
seqs_reads_boxplot<-
  compare_read_counts%>%
  mutate(Source=if_else(Source %in% c("Buccal","Saliva"),"Buccal/Saliva",Source))%>%
  mutate(Source=str_replace(Source,"Lung","Normal"))%>%
  mutate(Source=if_else(str_detect(Source,"Blank"),"Neg. Control",Source))%>%
  filter(!Source %in% c("Peritumoral")) %>%
  ggplot(aes(x=factor(Source,levels=c("Tumor","Normal","Neg. Control",
                                      "Blood","Buccal/Saliva")),
             y=log10(Reads+1),fill=dataset))+
  geom_boxplot(position=position_dodge2(preserve = "single"),
               outlier.alpha = 0,width=0.6)+
  coord_cartesian(ylim=c(1,6.5))+
  scale_fill_brewer(type = "qual")+
  labs(x=element_blank(),y="log10 Bacterial Reads")

## 1f
hmp_data <-
  read_delim("extra_data/HMP_data/hmps-braken-conf.1-collated.txt") %>%
  pivot_longer(-1, names_to = "sample_id", values_to = "Reads") %>%
  rename(tax_id = "4")
hmp_meta <-
  read_delim("extra_data/HMP_data/hmp_manifest_metadata_sourcetrack.tsv") %>%
  distinct() %>%
  left_join(read_delim("extra_data/HMP_data/hmp_manifest_final.tsv") %>% distinct()) %>%
  select(file_id, sample_body_site) %>%
  distinct()
hmp_readCounts <-
  hmp_data %>%
  left_join(hmp_meta, by = c("sample_id" = "file_id")) %>%
  filter(tax_id %in% all_bact_genera) %>%
  summarise(Reads = sum(Reads, na.rm = T),.by=c(sample_id, sample_body_site)) %>%
  mutate(sample_body_site = paste0("HMP ", sample_body_site)) %>%
  mutate(Source = "Re-analyzed")

# Numbers here for calculating reads/million come from the Gihawi paper
Gihawi_results <-
  read_delim("extra_data/BLCA_gihawi.txt") %>%
  pivot_longer(-1, names_to = "name", values_to = "Reads") %>%
  mutate(tissue = "Gihawi-PCAWG BLCA", TOTAL_READS = 205521556080/238) %>%
  rbind(read_delim("extra_data/BRCA_gihawi.txt") %>%
          pivot_longer(-1, names_to = "name", values_to = "Reads") %>%
          mutate_at("name", str_replace, "g_", "") %>%
          mutate(tissue = "Gihawi-PCAWG BRCA",TOTAL_READS = 324824097837/238)) %>%
  rbind(read_delim("extra_data/HNSC_gihawi.txt") %>%
          pivot_longer(-1, names_to = "name", values_to = "Reads") %>%
          mutate_at("name", str_replace, "g_", "") %>%
          mutate(tissue = "Gihawi-PCAWG HNSC",TOTAL_READS = 258961253944/334)) %>%
  # filter(!name %in% salter_list_nonHuman$Genus) %>%
  summarise(Reads = sum(Reads),.by=c(Sample, tissue,TOTAL_READS)) %>%
  mutate(Source = "Public",Reads_per_million=Reads/TOTAL_READS*1e6)%>%
  select(Sample,tissue,Source,TOTAL_READS,Reads,Reads_per_million)

wgs_reads_comparison <-
  kraken_raw %>%
  pivot_longer(-1, names_to = "Barcode", values_to = "Reads") %>%
  filter(
    Barcode %in% colnames(wgs_raw_combined) | str_detect(Barcode, "TCGA"),tax_id %in% all_bact_genera
  ) %>%
  left_join(annotations_withClinical) %>%
  #filter(!tax_id %in% salter_list_nonHuman$tax_id) %>%
  filter(Sample_Source == "Tumor") %>%
  mutate(Sample_Source = case_when(
    str_detect(Barcode, "TCGA") & Smoking == "Smoker" ~ "PCAWG Smoker LUAD",
    str_detect(Barcode, "TCGA") & Smoking == "Non-Smoker" ~ "PCAWG LCINS LUAD",
    T ~ "Sherlock WGS"
  )) %>%
  summarise(Reads = sum(Reads, na.rm = T),.by=c(Barcode, Sample_Source)) %>%
  mutate(Source = if_else(str_detect(Barcode, "TCGA"), "Re-analyzed", "Sherlock"))%>%
  left_join(wgs_read_depth)%>%
  drop_na(TOTAL_READS)%>%
  mutate(Reads_per_million=Reads/TOTAL_READS*1e6)%>%
  select(Barcode,Source,Sample_Source,TOTAL_READS,Reads,Reads_per_million)

rna_reads_comparison <-
  rna %>%
  pivot_longer(-1, names_to = "RNAseq_SampleID", values_to = "Reads") %>%
  filter(
    RNAseq_SampleID %in% colnames(rna_combatd_decontamd) | str_detect(RNAseq_SampleID, "TCGA"),
    tax_id %in% all_bact_genera
  ) %>%
  left_join(rna_annots[, c("RNAseq_SampleID", "Type", "SMOKING_HISTORY_txt")]) %>%
  #filter(!tax_id %in% salter_list_nonHuman$tax_id) %>%
  filter(Type == "Tumor") %>%
  mutate(Type = case_when(
    str_detect(RNAseq_SampleID, "TCGA") & SMOKING_HISTORY_txt == "smoker" ~ "TCGA LUAD/LUSC RNA",
    str_detect(RNAseq_SampleID, "TCGA") & SMOKING_HISTORY_txt == "never-smoker" ~ "TCGA LCINS RNA",
    SMOKING_HISTORY_txt == "never-smoker" ~ "Sherlock RNA",
    TRUE ~ "??"
  )) %>%
  summarise(Reads = sum(Reads, na.rm = T),.by=c(RNAseq_SampleID, Type)) %>%
  mutate(Source = if_else(str_detect(RNAseq_SampleID, "TCGA"), "Re-analyzed", "Sherlock"))%>%
  left_join(rna_total_reads)%>%
  mutate(Reads_per_million=Reads/as.numeric(TOTAL_READS)*1e6)%>%
  select(RNAseq_SampleID,Source,Type,TOTAL_READS,Reads,Reads_per_million)

multiOrgan_comparison_fig <-
  Gihawi_results %>%
  mutate(Source = "Public",dataset="WGS") %>%
  #rbind(hmp_readCounts %>% rename(Sample = "sample_id", tissue = "sample_body_site")) %>%
  rbind(wgs_reads_comparison %>% rename(Sample = "Barcode", tissue = "Sample_Source")%>%mutate(dataset="WGS")) %>%
  rbind(rna_reads_comparison %>% rename(Sample = "RNAseq_SampleID",tissue="Type")%>%mutate(dataset="RNA-seq")) %>%
  filter(!tissue == "??") %>%
  ggplot(aes(
    x = fct_reorder(tissue, .x = Reads_per_million, .fun = median, .desc = T),
    y = log10(Reads_per_million), fill = Source
  )) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.075), alpha = 0.075) +
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45)) +
  labs(x = NULL,y="log10 Bacterial Reads Per Million") +
  scale_fill_brewer(type="qual",palette = "Set2")+
  facet_grid2(~dataset,scales = "free_x",
              strip = strip_themed(background_x = list(element_rect(fill=RNA_color),
                                                       element_rect(fill=WGS_color))),
              space ="free_x")

fig1_legend <- get_legend(kingdom_barplot)

plot_grid(
  plot_grid(sankey_plot, workflow, ncol = 1, rel_heights = c(.66, .33),labels=c("a","b")),
  plot_grid(kingdom_barplot,
            plot_grid(bacterial_proportions_boxplot + 
                        guides(fill = "none") + theme(axis.title.y = element_text(size = 10)),
                      seqs_reads_boxplot + guides(fill = "none"),
                      fig1_legend,labels=c("d","e",""),rel_widths = c(1, 1,.5),nrow=1),
            multiOrgan_comparison_fig,labels = c("c","","f"),
            nrow=3,ncol=1,align="hv",rel_heights = c(0.4,0.3,0.4)),
  labels = c("a","",""), rel_widths = c(0.45, 0.55))

ggsave("plots/current_plots/fig1_v1.pdf",
       device = cairo_pdf, width = 11 / 1.1, height = 10 / 1.1, bg = "white"
)
ggsave("plots/current_plots/fig1_v1.png",
       device = png, width = 11 / 1.1, height = 10 / 1.1, bg = "white"
)