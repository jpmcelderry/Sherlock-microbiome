####### FIGURE 2

######## Big phylum-level compositions plot
graph_phyla_composition_barplot2 <-
  function(pivoted_relabund_dataset, ordering_dataset = pivoted_relabund_dataset, level = "phylum",
           ordering_taxa = "Other", otus = 5, relabundance_cutoff = 0,
           reads_cutoff = 1, legend_size = 8, decreasing = F) {
    if (otus > nlevels(as.factor(pivoted_relabund_dataset$name))) {
      otus <- nlevels(as.factor(pivoted_relabund_dataset$name))
    }
    tmp <-
      pivoted_relabund_dataset %>%
      filter(type == level) %>%
      filter(Relabund >= relabundance_cutoff, Reads >= reads_cutoff) %>%
      mutate(lumped = fct_lump(name, n = otus))
    ordering <-
      tmp %>%
      select(Barcode, lumped, Relabund) %>%
      complete(Barcode, lumped, fill = list(Relabund = 0)) %>%
      filter(lumped == ordering_taxa) %>%
      group_by(Barcode) %>%
      summarise(Relabund = sum(Relabund)) %>%
      arrange((1 - 2 * decreasing) * Relabund) %>%
      pull(Barcode)
    tmp_plot <-
      tmp %>%
      ggplot(aes(
        y = Relabund,
        x = factor(Barcode, levels = ordering),
        fill = lumped
      )) +
      geom_bar(width = 1, position = "fill", stat = "identity") +
      scale_fill_manual(values = SBScolor) +
      scale_x_discrete(labels = NULL) +
      theme_bw(base_family = "Roboto Condensed") +
      theme(
        legend.text = element_text(size = legend_size),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      ) +
      labs(fill = level)
    return(tmp_plot)
  }
wgs_compOverview_in <-
  wgs_raw_combined %>%
  as.data.frame() %>%
  filter(tax_id %in% all_bact_phyla) %>%
  right_join(kraken_taxonomy, .) %>%
  pivot_longer(-c(1:4), names_to = "Barcode", values_to = "Reads") %>%
  mutate(Relabund = Reads / sum(Reads),.by=Barcode) %>%
  # select(-Barcode)%>%
  left_join(wgs_full_annotations[, c("Barcode", "Sample_Source", "sherlock_pid")],
            by = c("Barcode")
  ) %>%
  filter(Sample_Source %in% c("Tumor", "Normal", "Blood")) %>%
  mutate(experiment = paste0("WGS (n=", nlevels(as.factor(sherlock_pid)), " subjects)"), .after = "tax_id")
s16_compOverview_in <-
  s16_scrubbed %>%
  as.data.frame() %>%
  # mutate(experiment=paste0("16s (n=",ncol(.)-1,")"),.after="tax_id")%>%
  filter(tax_id %in% all_bact_phyla) %>%
  right_join(kraken_taxonomy, .) %>%
  pivot_longer(-c(1:4), names_to = "Barcode", values_to = "Reads") %>%
  mutate(Relabund = Reads / sum(Reads),.by=Barcode) %>%
  # select(-Barcode)%>%
  left_join(s16_metadata[, c("SampleID", "Tumor-NormalStatus", "AdditionalAttributes")],
            by = c("Barcode" = "SampleID")
  ) %>%
  filter(`Tumor-NormalStatus` %in% c("Tumor", "Normal")) %>%
  mutate(experiment = paste0("16S (n=", nlevels(as.factor(AdditionalAttributes)), " subjects)"), .after = "tax_id")
RNAseq_compOverview_in <-
  rna_rawCounts_decontamd %>%
  as.data.frame() %>%
  filter(tax_id %in% all_bact_phyla) %>%
  right_join(kraken_taxonomy, .) %>%
  pivot_longer(-c(1:4), names_to = "Barcode", values_to = "Reads") %>%
  mutate(Relabund = Reads / sum(Reads),.by=Barcode) %>%
  left_join(rna_annots[, c("RNAseq_SampleID", "Type", "Subject ID")],
            by = c("Barcode" = "RNAseq_SampleID")
  ) %>%
  filter(Type %in% c("Tumor", "Normal")) %>%
  mutate(experiment = paste0("RNA-seq (n=", nlevels(as.factor(`Subject ID`)), " subjects)"), .after = "tax_id")
composition_plot_grid <-
  RNAseq_compOverview_in %>%
  rbind(s16_compOverview_in %>% rename(Type = "Tumor-NormalStatus", `Subject ID` = "AdditionalAttributes")) %>%
  rbind(wgs_compOverview_in %>% rename(Type = "Sample_Source", `Subject ID` = "sherlock_pid")) %>%
  mutate(Type = factor(Type, levels = c("Normal", "Tumor", "Blood"))) %>%
  mutate(name=case_when(name=="Pseudomonadota"~"Proteobacteria",
                        name=="Actinomycetota"~"Actinobacteria",
                        name=="Bacillota"~"Firmicutes",
                        name=="Bacteroidota"~"Bacteroidetes",
                        name=="Deinococcota"~"Deinoccocus-Thermus",
                        TRUE~name
  ))%>%
  graph_phyla_composition_barplot2(level = "phylum", ordering_taxa = "Proteobacteria", otus = 6) +
  facet_grid2(Type ~ experiment,
              scales = "free",
              independent = "x",
              render_empty = F,
              strip = strip_themed(background_x = list(element_rect(fill=s16_color),
                                                       element_rect(fill=RNA_color),
                                                       element_rect(fill=WGS_color)))
  ) +
  facetted_pos_scales(y = list(Type == "Blood" ~ scale_y_continuous(breaks = NULL))) +
  theme(
    legend.position = c(0, 0),
    legend.justification = c(0, -.5),
    legend.title.position = "top",
    legend.byrow = T,
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(12, units = "pt")
  ) +
  guides(fill = guide_legend(nrow = 2)) +
  # guides(color = guide_legend(override.aes = list(size = 1))) +
  scale_fill_manual(values = SBScolor[c(1:4,7:12)]) +
  labs(x = NULL, y = "Relative Abundance", fill = "Phylum")
composition_plot_grid


##### Composition summary barplots
wgs_lung_samples<-
  annotations_withClinical%>%
  filter(Sample_Source %in% c("Lung","Tumor"))%>%
  pull(Barcode)

# Phylum level
relabund_byExperiment2<-
  rna_rawCounts_decontamd%>% # Join data
  select(tax_id,any_of(samples_crossPlatform$Barcode))%>%
  mutate(experiment="RNA-seq",.after="tax_id")%>%
  full_join(wgs_raw_combined%>%
              mutate(experiment="WGS",.after="tax_id"))%>%
  full_join(s16_scrubbed%>%
              mutate(experiment="16S",.after="tax_id"))%>%
  right_join(unified_taxonomy,.)%>%
  filter(tax_id %in% all_bact_phyla)%>%
  pivot_longer(-c(1:5),names_to = "Sample",values_to = "Reads")%>% # Pivot data, minus taxonomy info
  replace_na(list(Reads=0))%>%
  mutate(Relabund=Reads/sum(Reads),.by=Sample)%>% # Transform to relative abundance
  rename(Barcode="Sample")%>%
  filter(Barcode %in% samples_crossPlatform$Barcode)%>% # Filter to only cross-platform samples
  left_join(all_sample_TNstatus[,c("Barcode","Source")])%>% # Join Tumor-Normal info
  drop_na()%>%
  filter(!(Source=="Blood"))%>%
  summarise(Relabund=mean(Relabund),.by=c(experiment,Source,name))%>% # Determine Platform/Tissue Type mean relAbund
  mutate(rank=order(order(Relabund, decreasing=TRUE)),.by=c(experiment,Source))%>% # Determine relAbund ranks per sample group
  mutate(min_rank=min(rank),.by=name) # Calculate the lowest rank per bacteria among all platforms/tissues 
phylum_lvl<-
  relabund_byExperiment2%>%
  mutate(name=if_else(min_rank>5,"Other",name))%>% # Lump any bacteria that's never top 5 abundant
  summarise(Relabund=sum(Relabund),.by=c(name,Source,experiment))%>% # Sum the "Other"s into 1 datapoint
  mutate(name=case_when(name=="Pseudomonadota"~"Proteobacteria",
                        name=="Actinomycetota"~"Actinobacteria",
                        name=="Bacillota"~"Firmicutes",
                        name=="Bacteroidota"~"Bacteroidetes",
                        name=="Deinococcota"~"Deinoccocus-\nThermus",
                        TRUE~name
  ))%>%
  mutate(name=fct_relevel(name,"Other",after=Inf))%>% # Place Other as last category
  mutate(nsamples=n())%>%
  ggplot(aes(x=Source,y=Relabund,fill=name))+
  geom_bar(position = "fill",stat="identity")+
  scale_fill_manual("Phylum",values = SBScolor[c(1:4,7:12)])+
  labs(y="Mean Relative Abundance",x="n=798")+
  theme(strip.placement = "outside",
        panel.spacing=unit(0,"cm"),
        axis.text.x = element_text(vjust=1,hjust=1,angle=45))+
  guides(color=guide_legend(override.aes = list(size = 10)))+
  ggh4x::facet_grid2(~experiment, space="free", scales="free",
                     strip = strip_themed(background_x = list(element_rect(fill=s16_color),
                                                              element_rect(fill=RNA_color),
                                                              element_rect(fill=WGS_color))))
phylum_lvl

# Genus level summary barplot
relabund_byExperiment<-
  rna_rawCounts_decontamd%>% # Join data
  mutate(experiment="RNA-seq",.after="tax_id")%>%
  full_join(wgs_raw_combined%>%
              mutate(experiment="WGS",.after="tax_id"))%>%
  full_join(s16_scrubbed%>%
              mutate(experiment="16S",.after="tax_id"))%>%
  right_join(unified_taxonomy,.)%>%
  filter(tax_id %in% all_bact_genera)%>%
  pivot_longer(-c(1:5),names_to = "Sample",values_to = "Reads")%>% # Pivot data, minus taxonomy info
  replace_na(list(Reads=0))%>%
  mutate(Relabund=Reads/sum(Reads),.by=Sample)%>% # Transform to relative abundance
  rename(Barcode="Sample")%>%
  filter(Barcode %in% samples_crossPlatform$Barcode)%>% # Filter to only cross-platform samples
  left_join(all_sample_TNstatus[,c("Barcode","Source")])%>% # Join Tumor-Normal info
  drop_na()%>%
  filter(!(Source=="Blood"))%>%
  summarise(Relabund=mean(Relabund),.by=c(experiment,Source,name))%>% # Determine Platform/Tissue Type mean relAbund
  mutate(rank=order(order(Relabund, decreasing=TRUE)),.by=c(experiment,Source))%>% # Determine relAbund ranks per sample group
  mutate(min_rank=min(rank),.by=name) # Calculate the lowest rank per bacteria among all platforms/tissues 
bottom_lvl<-
  relabund_byExperiment%>%
  mutate(name=if_else(min_rank>12,"Other",name))%>% # Lump any bacteria that's never top 15 abundant
  summarise(Relabund=sum(Relabund),.by=c(name,Source,experiment))%>% # Sum the "Other"s into 1 datapoint 
  mutate(name=fct_relevel(name,"Other",after=Inf))%>%# Place Other as last category
  ggplot(aes(x=Source,y=Relabund,fill=name))+
  geom_bar(position = "fill",stat="identity")+
  scale_fill_manual("Genus",values=unname(SBScolor))+
  labs(y="Mean Relative Abundance",x="n=798")+
  facet_grid2(~experiment, space="free_x", scales="free_x",
              strip = strip_themed(background_x = list(element_rect(fill=s16_color),
                                                       element_rect(fill=RNA_color),
                                                       element_rect(fill=WGS_color))))+
  theme(strip.placement = "outside",
        panel.spacing=unit(0,"cm"),
        axis.text.x = element_text(vjust=1,hjust=1,angle=45))
bottom_lvl

# genus_phylum_info<-
#   tax_table%>%
#   rownames_to_column("tax_id")%>%
#   drop_na(genus)%>%
#   select(genus=tax_id,phylum)%>%
#   left_join(kraken_taxonomy,by=c("phylum"="tax_id"))%>%
#   select(genus,phylum,phylum_name=name)


##### species level?
relabund_byExperiment3<-
  rna%>% # Join data
  mutate(experiment="RNA-seq",.after="tax_id")%>%
  full_join(read_delim("/Volumes/Sherlock_Lung/JohnMce/Kraken_WGS3/kraken0.1-collated.txt")%>%
              rename(tax_id=`4`)%>%
              mutate(tax_id=as.character(tax_id))%>%
              mutate(experiment="WGS",.after="tax_id"))%>%
  right_join(unified_taxonomy,.)%>%
  filter(type=="species",str_detect(taxonomy,"Bacteria"))%>%
  mutate(genus=str_split(name," ")%>%sapply(`[`,1),.after = name)%>%
  filter(!genus %in% c(salter_list_nonHuman$Genus,"Cutibacterium"))%>%
  pivot_longer(-c(1:6),names_to = "Sample",values_to = "Reads")%>% # Pivot data, minus taxonomy info
  replace_na(list(Reads=0))%>%
  mutate(Relabund=Reads/sum(Reads),.by=Sample)%>% # Transform to relative abundance
  rename(Barcode="Sample")%>%filter(Barcode %in% samples_crossPlatform$Barcode)%>% # Filter to only cross-platform samples
  left_join(all_sample_TNstatus[,c("Barcode","Source")])%>% # Join Tumor-Normal info
  drop_na()%>%
  filter(!(Source=="Blood"))%>%
  summarise(Relabund=sum(Relabund)/n(),.by=c(experiment,Source,name,tax_id))%>% # Determine Platform/Tissue Type mean relAbund
  mutate(rank=order(order(Relabund, decreasing=TRUE)),.by=c(experiment,Source))%>% # Determine relAbund ranks per sample group
  mutate(min_rank=min(rank),.by=name) # Calculate the lowest rank per bacteria among all platforms/tissues 
species_lvl<-
  relabund_byExperiment3%>%
  mutate(name=if_else(min_rank>20,"Other",name))%>% # Lump all bacteria that's never 'top 25' abundant
  summarise(Relabund=sum(Relabund),.by=c(name,Source,experiment))%>% # Sum the "Other"s into 1 datapoint 
  mutate(name=fct_relevel(name,"Other",after=Inf))%>%# Place Other as last category
  ggplot(aes(x=Source,y=Relabund,fill=name))+
  geom_bar(position = "fill",stat="identity")+
  scale_fill_manual("Genus",values=unname(SBScolor))+
  labs(y="Mean Relative Abundance",x="n=798")+
  guides(fill=guide_legend(ncol=2))+
  facet_grid2(~experiment, space="free_x", scales="free_x",
              strip = strip_themed(background_x = list(element_rect(fill=s16_color),
                                                       element_rect(fill=RNA_color),
                                                       element_rect(fill=WGS_color))))+
  theme(strip.placement = "outside",
        panel.spacing=unit(0,"cm"),
        axis.text.x = element_text(vjust=1,hjust=1,angle=45))
species_lvl

# Fig 2d
big_rarecurve<-
  rna_rawCounts_decontamd%>%
  full_join(wgs_raw_combined)%>%
  full_join(s16_scrubbed)%>%
  filter(tax_id %in% all_bact_genera)%>%
  replace(is.na(.),0)%>%
  column_to_rownames("tax_id")%>%
  t()%>%
  rarecurve(step=100,MARGIN=2,tidy = TRUE)

rarefaction_plot<-
  big_rarecurve%>%
  rename(Barcode="Site")%>%
  left_join(all_sample_TNstatus)%>%
  mutate(dataset=str_replace(dataset,pattern = "RNA",replacement = "RNA-seq"))%>%
  filter(Sample<=2001,!Source=="Peritumoral")%>% ### Cap the sampling to 2,500 reads
  ggplot(aes(x=Sample,y=Species,color=dataset))+
  #geom_line(aes(group=Barcode),alpha=0.05)+
  geom_smooth(method = "loess")+
  labs(y="Genus richness")+
  facet_wrap(~Source,ncol=2)+
  theme(legend.position = c(.9, 0), 
        legend.justification = c(1,-.1),
        axis.text.x = element_text(size=7))+
  scale_color_brewer(type = "qual")+
  labs(x="Rarefaction depth")

# plot Figure 2
plot_grid(composition_plot_grid,
          plot_grid(phylum_lvl +
                      theme(
                        legend.text = element_text(size = 8),
                        legend.title = element_text(size = 10),
                        legend.key.size = unit(8, units = "pt"),
                      ),
                    bottom_lvl +
                      theme(
                        legend.text = element_text(size = 8),
                        legend.title = element_text(size = 10),
                        legend.key.size = unit(8, units = "pt")
                      ) +
                      guides(fill = guide_legend(ncol = 1)),
                    rarefaction_plot,
                    labels = c("b", "c", "d"),
                    nrow = 1, align = "h",
                    rel_widths = c(1.3, 1.4, 1)),
          labels = c("a",""),
          nrow=2,
          rel_heights = c(1.5,1))
ggsave("plots/current_plots/fig2_v4.pdf", device = cairo_pdf, height = 10, width = 11)
ggsave("plots/current_plots/fig2_v4.png", device = png, height = 10, width = 11)

