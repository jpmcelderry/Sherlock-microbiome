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
  right_join(unified_taxonomy, .) %>%
  filter(type %in% c("phylum","genus"),str_detect(taxonomy,"Bacteria")) %>%
  pivot_longer(-c(1:4), names_to = "Barcode", values_to = "Reads") %>%
  mutate(Relabund = Reads / sum(Reads),.by=c(Barcode,type)) %>%
  left_join(wgs_full_annotations[, c("Barcode", "Sample_Source", "sherlock_pid")],
            by = c("Barcode")
  ) %>%
  filter(Sample_Source %in% c("Tumor", "Normal", "Blood")) %>%
  mutate(experiment = paste0("WGS (n=", nlevels(as.factor(sherlock_pid)), " subjects)"), .after = "tax_id")
s16_compOverview_in <-
  s16_scrubbed %>%
  right_join(unified_taxonomy, .) %>%
  filter(type %in% c("phylum","genus"),str_detect(taxonomy,"Bacteria")) %>%
  pivot_longer(-c(1:4), names_to = "Barcode", values_to = "Reads") %>%
  mutate(Relabund = Reads / sum(Reads),.by=c(Barcode,type)) %>%
  left_join(s16_metadata[, c("SampleID", "Tumor-NormalStatus", "AdditionalAttributes")],
            by = c("Barcode" = "SampleID")
  ) %>%
  filter(`Tumor-NormalStatus` %in% c("Tumor", "Normal")) %>%
  mutate(experiment = paste0("16S (n=", nlevels(as.factor(AdditionalAttributes)), " subjects)"), .after = "tax_id")
RNAseq_compOverview_in <-
  rna_rawCounts_decontamd %>%
  right_join(unified_taxonomy, .) %>%
  filter(type %in% c("phylum","genus"),str_detect(taxonomy,"Bacteria")) %>%
  pivot_longer(-c(1:4), names_to = "Barcode", values_to = "Reads") %>%
  mutate(Relabund = Reads / sum(Reads),.by=c(Barcode,type)) %>%
  left_join(rna_annots[, c("RNAseq_SampleID", "Type", "Subject ID")],
            by = c("Barcode" = "RNAseq_SampleID")
  ) %>%
  filter(Type %in% c("Tumor", "Normal")) %>%
  mutate(experiment = paste0("RNA-seq (n=", nlevels(as.factor(`Subject ID`)), " subjects)"), .after = "tax_id")

composition_plot_grid <-
  RNAseq_compOverview_in %>%
  filter(type=="phylum")%>%
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
#rm(wgs_compOverview_in,s16_compOverview_in,RNAseq_compOverview_in)

##### Composition summary barplots
# Compile samples shared across platforms
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
  filter(WGS %in% colnames(wgs_decontamd)|WGS %in% colnames(NG232_decontamd),
         `16s` %in% colnames(s16_scrubbed),
         RNA %in% colnames(rna_combatd_decontamd))%>%
  pivot_longer(c("WGS","16s","RNA"),names_to = "platform",values_to = "Barcode")

#
relabund_byExperiment<-
  # Join data
  rna_rawCounts_decontamd%>% 
  select(tax_id,any_of(samples_crossPlatform$Barcode))%>%
  mutate(experiment="RNA-seq",.after="tax_id")%>%
  full_join(wgs_raw_combined%>%
              mutate(experiment="WGS",.after="tax_id"))%>%
  full_join(s16_scrubbed%>%
              mutate(experiment="16S",.after="tax_id"))%>%
  right_join(unified_taxonomy,.)%>%
  filter(type %in% c("phylum","genus"),str_detect(taxonomy,"Bacteria"))%>%
  # Transform to relative abundance
  pivot_longer(-c(1:5),names_to = "Sample",values_to = "Reads")%>%
  replace_na(list(Reads=0))%>%
  mutate(Relabund=Reads/sum(Reads),.by=c(Sample,type))%>% 
  rename(Barcode="Sample")%>%
  # Filter to only cross-platform samples, add TN information
  filter(Barcode %in% samples_crossPlatform$Barcode)%>% 
  left_join(all_sample_TNstatus[,c("Barcode","Source")])%>%
  drop_na()%>%
  # add info for lumping least abundant bacteria downstream
  summarise(Relabund=mean(Relabund),n=n()/3,.by=c(experiment,Source,name,type))%>% 
  mutate(rank=order(order(Relabund, decreasing=TRUE)),.by=c(experiment,Source,type))%>% 
  mutate(highest_rank=min(rank),.by=c(name,type))

genus_level_stackedBar<-
  relabund_byExperiment%>%
  filter(type=="genus")%>%
  # Lump bacteria into 'Other' group and calculate it's total abundance
  mutate(name=if_else(highest_rank>12,"Other",name))%>%
  summarise(Relabund=sum(Relabund),.by=c(name,Source,experiment,n))%>%
  mutate(name=fct_relevel(name,"Other",after=Inf))%>%
  # plot
  ggplot(aes(x=Source,y=Relabund,fill=name))+
  geom_bar(position = "fill",stat="identity")+
  scale_fill_manual("Genus",values=unname(SBScolor))+
  labs(y="Mean Relative Abundance",x="n=278 Normal,522 Tumor")+
  facet_grid2(~experiment, space="free_x", scales="free_x",
              strip = strip_themed(background_x = list(element_rect(fill=s16_color),
                                                       element_rect(fill=RNA_color),
                                                       element_rect(fill=WGS_color))))+
  theme(strip.placement = "outside",
        panel.spacing=unit(0,"cm"),
        axis.text.x = element_text(vjust=1,hjust=1,angle=45))

phylum_level_stackedBar<-
  relabund_byExperiment%>%
  filter(type=="phylum")%>%
  # Lump bacteria into 'Other' group and calculate it's total abundance
  mutate(name=if_else(highest_rank>5,"Other",name))%>%
  summarise(Relabund=sum(Relabund),.by=c(name,Source,experiment,n))%>%
  mutate(name=case_when(name=="Pseudomonadota"~"Proteobacteria",
                        name=="Actinomycetota"~"Actinobacteria",
                        name=="Bacillota"~"Firmicutes",
                        name=="Bacteroidota"~"Bacteroidetes",
                        name=="Deinococcota"~"Deinococcus-\nThermus",
                        TRUE~name
  ))%>%
  mutate(name=fct_relevel(name,"Other",after=Inf))%>% # Place Other as last category
  mutate(nsamples=n())%>%
  # plot
  ggplot(aes(x=Source,y=Relabund,fill=name))+
  geom_bar(position = "fill",stat="identity")+
  scale_fill_manual("Phylum",values = SBScolor[c(1:4,7:12)])+
  labs(y="Mean Relative Abundance",x="n=278 Normal,522 Tumor")+
  theme(strip.placement = "outside",
        panel.spacing=unit(0,"cm"),
        axis.text.x = element_text(vjust=1,hjust=1,angle=45))+
  guides(color=guide_legend(override.aes = list(size = 10)))+
  ggh4x::facet_grid2(~experiment, space="free", scales="free",
                     strip = strip_themed(background_x = list(element_rect(fill=s16_color),
                                                              element_rect(fill=RNA_color),
                                                              element_rect(fill=WGS_color))))

# Fig 2d
big_rarecurve<-
  # calculate rarefaction curve
  rna_rawCounts_decontamd%>%
  full_join(wgs_raw_combined)%>%
  full_join(s16_scrubbed)%>%
  filter(tax_id %in% all_bact_genera)%>%
  replace(is.na(.),0)%>%
  column_to_rownames("tax_id")%>%
  t()%>%
  rarecurve(step=100,MARGIN=2,tidy = TRUE)

rarefaction_plot<-
  ### Cap the sampling to 2,000 reads, merge metadata
  big_rarecurve%>%
  filter(Sample<=2001)%>%
  rename(Barcode="Site")%>%
  left_join(all_sample_TNstatus)%>%
  filter(!Source=="Peritumoral")%>%
  mutate(dataset=str_replace(dataset,pattern = "RNA",replacement = "RNA-seq"))%>%
  # plot
  ggplot(aes(x=Sample,y=Species,color=dataset))+
  geom_smooth(method = "loess")+
  labs(y="Richness")+
  facet_wrap(~Source,ncol=2)+
  theme(legend.position = c(.9, 0), 
        legend.justification = c(1,-.1),
        axis.text.x = element_text(size=7))+
  scale_color_brewer(type = "qual")+
  labs(x="Rarefaction depth")

# plot Figure 2
plot_grid(composition_plot_grid,
          plot_grid(phylum_level_stackedBar +
                      theme(
                        legend.text = element_text(size = 8),
                        legend.title = element_text(size = 10),
                        legend.key.size = unit(8, units = "pt"),
                      ),
                    genus_level_stackedBar +
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

