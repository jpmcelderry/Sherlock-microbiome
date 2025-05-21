# Extended Data Fig 1 - Genus level composition barplot
composition_plot_grid_genera <-
  # MERGE DATASETS
  RNAseq_compOverview_in %>%
  filter(type=="genus")%>%
  rbind(s16_compOverview_in %>% rename(Type = "Tumor-NormalStatus", `Subject ID` = "AdditionalAttributes")) %>%
  rbind(wgs_compOverview_in %>% rename(Type = "Sample_Source", `Subject ID` = "sherlock_pid")) %>%
  mutate(Type = factor(Type, levels = c("Normal", "Tumor", "Blood"))) %>%
  graph_phyla_composition_barplot2(level = "genus", ordering_taxa = "Other", otus = 20) +
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
  guides(fill = guide_legend(nrow = 6)) +
  # guides(color = guide_legend(override.aes = list(size = 1))) +
  scale_fill_manual(values = SBScolor) +
  labs(x = NULL, y = "Relative Abundance", fill = "Genus")
composition_plot_grid_genera

ggsave("plots/current_plots/composition_panels_genera.pdf",device = cairo_pdf,height=8,width=9)
ggsave("plots/current_plots/composition_panels_genera.png",device = png,height=8,width=9)

plot_grid(richness_v_immuneCells_plot+guides(color=FALSE,size=FALSE),
                    diversity_v_immuneCells_plot+guides(color=FALSE,size=FALSE),
                    get_legend(diversity_v_immuneCells_plot+guides(color=guide_legend(override.aes = list(alpha=1)))),
                    s16_immune_richness_plot+guides(color=FALSE,size=FALSE),
                    s16_immune_diversity_plot+guides(color=FALSE,size=FALSE),
                    get_legend(s16_immune_diversity_plot+guides(color=guide_legend(override.aes = list(alpha=1)))),
                    labels=c("a","","","b","",""),nrow=2,rel_widths = c(.5,.5,.2,.5,.5,.2))
ggsave("plots/current_plots/immune_cells_supplement.pdf",device=cairo_pdf,width=9,height=8,bg="white")
ggsave("plots/current_plots/immune_cells_supplement.png",device=png,width=9,height=8,bg="white")

# Extended Data Fig 3 - Genomics associations
source("genomics_associations.R")

######## batch correction, supplemental figure 1
WGS_batchCorrect_figures<-
  plot_grid(
    plot_grid(
      WGS_batchBefore_plot+theme(legend.position = "none"),
      plot_grid(
        WGS_batchAfter_plot+theme(legend.position = "none"),
        NG232_batchAfter_plot+theme(legend.position = "none"),
        labels=c("b","c"),
        rel_widths = c(1,1)
      ),
      labels = c("a"),rel_widths = c(1,1),nrow=2
    ),
    get_legend(WGS_batchBefore_plot),
    nrow=1,ncol=2,rel_widths=c(1.5,.5)
  )
RNA_batchCorrect_figures<-
  plot_grid(RNA_batchBefore_plot,
            RNA_batchAfter_plot+theme(legend.position = "none"),
            get_legend(RNA_batchAfter_plot),
            labels = c("d","e",""),rel_widths = c(1,1,.5),nrow=1)
s16_batchCorrect_figures<-
  plot_grid(s16_batchAfter_plot+theme(legend.position = "none"),
            get_legend(s16_batchAfter_plot),
            NULL,labels=c("f",""),rel_widths = c(1,.5,1),
            nrow=1)

plot_grid(WGS_batchCorrect_figures,
          RNA_batchCorrect_figures,
          s16_batchCorrect_figures,
          ncol=1,rel_heights = c(2,1,1))
ggsave("plots/current_plots/batch_correct.pdf",device=cairo_pdf,width=10,height=11,bg="white")
ggsave("plots/current_plots/batch_correct.png",device=png,width=10,height=11,bg="white")

# Supplementary Fig 2 - Raw data description
source("raw_composition_supp.R")

# Supplementary Fig 3 - ICC metrics
source("ICC_metrics.R")
# Relabunance correlations
fs_ICCs_top<-
  plot_grid(WGS_16s_phylum_correlation,
            RNA_16s_phylum_correlation,
            RNA_WGS_phylum_correlation,
            labels=c("a","c","e"),
            ncol=1,rel_heights = c(0.3,0.45,0.3))
# Diversity correlations
fs_ICCs_bottom<-
  plot_grid(wgs_s16_diversity,
            rna_s16_diversity,
            rna_wgs_diversity,
            labels=c("b","d","f"),
            ncol=1)
plot_grid(fs_ICCs_top,fs_ICCs_bottom,rel_widths = c(0.6,0.4))
ggsave("plots/current_plots/ICC_correlations.pdf",device=cairo_pdf,width=9,height=8)
ggsave("plots/current_plots/ICC_correlations.png",width=9,height=8)

# Supplementary Figure 4 - species differential abundance, supplementary diversity info
source("species_analysis.R")

# betadiversity plot
betaresults <-
  # bind tables
  tidy(big_model_16s)%>% 
    select(term, p.value, R2) %>%
    mutate(fdr = p.adjust(p.value, method = "fdr"),Experiment = "16S") %>%
  rbind(
    tidy(big_model_rna2) %>% 
      select(term, p.value, R2) %>% 
      mutate(fdr = p.adjust(p.value, method = "fdr"),Experiment = "RNA-seq")) %>%
  rbind(
    wgs_betadiver_meta_results %>% 
      mutate(Experiment = "WGS") %>% 
      rename(fdr="meta.fdr",p.value="meta.p.value")) %>%
  filter(!term %in% c("Residual", "Total")) %>%
  #mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  # Pretty up labels
  mutate(term = case_when(
    term %in% c("site_study", "Site_REF", "Study") ~ "Study Site",
    term == "Sample_Source" ~ "Tumor/Normal",
    term == "HISTOLOGY_COMPOSITE" ~ "Histology",
    term == "SEX_DERIVED" ~ "Sex",
    term == "STAGE_simple" ~ "Stage",
    term == "METASTASIS" ~ "Metastasis",
    term == "GRADE_DIFFERENTIATION" ~ "Grade/Differentiation",
    term == "ANY_PREVIOUS_LUNG_DISEASE" ~ "Any Previous Lung Disease",
    term == "PASSIVE_SMOKE" ~ "Passive Smoke",
    term == "ANCESTRY_DERIVED" ~ "Ancestry",
    term == "RECURRENCE" ~ "Recurrence",
    term == "VITAL_STATUS" ~ "Vital Status",
    term == "AGE_AT_DIAGNOSIS" ~ "Age at Diagnosis",
    TRUE ~ term
  )) %>%
  # plot
  ggplot(aes(
    x = term, y = R2, 
    fill = factor(Experiment, levels = c("16S", "RNA-seq","WGS")),
    group = factor(Experiment, levels = c("16S", "RNA-seq","WGS"))
  )) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  geom_text(
    aes(label = case_when(fdr < 0.0001 ~ "****",
                          fdr < 0.001 ~ "***",
                          fdr < 0.01 ~ "**",
                          fdr < 0.05 ~ "*",
                          .default = ""
    ), y = (R2 + .01 * R2)),
    position = position_dodge(width = .6), size = 12 / .pt
  ) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) +
  scale_fill_brewer(type = "qual") +
  labs(fill = "Experiment", y = "Proportion of Variance Explained", x = NULL)

plot_grid(
  plot_grid(aldex_species,
            rna_richness_country,
            labels = c("a","b")),
  plot_grid(hospital_v_aDiversity,
            betaresults+theme(legend.position = "top"),
            labels = c("c","d"),nrow=1),
  ncol=1,align="hv"
)
ggsave(filename = "plots/current_plots/suppFigure_diveristyClinical.pdf",device=cairo_pdf,width=8,height=8,bg = "white")
ggsave(filename = "plots/current_plots/suppFigure_diveristyClinical.png",device=png,width=8,height=8,bg = "white")

# Supplementary Fig 4 - blood diversity results
source("WGS_blood_analyses.R")

# Supplementary Figure 5 - see Supp Figure 1 code

# Supplementary Figure 6 - richness/diversity survival
# calculated in *_survival.R scripts
plot_grid(survivalRichness+
            theme(legend.text = element_text(size=10),legend.title = element_text(size=10)),
          survivalDiversity+
            theme(legend.text = element_text(size=10),legend.title = element_text(size=10)),
          survivalRichness_16S+
            theme(legend.text = element_text(size=10),legend.title = element_text(size=10)),
          survivalDiversity_16s+
            theme(legend.text = element_text(size=10),legend.title = element_text(size=10)),
          survivalRichness_wgs+
            theme(legend.text = element_text(size=10),legend.title = element_text(size=10)),
          survivalDiversity_wgs+
            theme(legend.text = element_text(size=10),legend.title = element_text(size=10)),
          Blood_survivalRichness_wgs+
            theme(legend.text = element_text(size=10),legend.title = element_text(size=10)),
          Blood_survivalDiversity_wgs+
            theme(legend.text = element_text(size=10),legend.title = element_text(size=10)),
          labels = c("a","b","c","d","e","f","g","h"),
          ncol=2)
ggsave(filename = "plots/current_plots/survival_diversity.pdf",device = cairo_pdf,height=9,width=8,bg="white")
ggsave(filename = "plots/current_plots/survival_diversity.png",device = png,height=9,width=8,bg="white")

# Supplementary Figure 7 - LUAD richness/diversity survival calculated with figure 4

#### Supplementary Figure 8 - CONFIDENCE SCORE DATA
testing_confidence<-
  read_delim('kraken_out/simulation-50000-all.txt',col_names = F)
testing_confidenceSpec<-
  testing_confidence%>%
  mutate(X3=str_split(X3," \\(taxid "))%>%
  unnest_wider(X3,names_sep = "classification")%>%
  mutate(X3classification2=str_replace(X3classification2,"\\)",""))%>%
  left_join(kraken_taxonomy%>%mutate(tax_id=as.character(tax_id)),by=c("X3classification2"="tax_id"))%>%
  filter(type=="species")
true_positives<-c("Escherichia coli","Pseudomonas aeruginosa",
                  "Klebsiella pneumoniae","Prevotella melaninogenica",
                  "Rothia mucilaginosa","Haemophilus parainfluenzae",
                  "Cutibacterium acnes","Moraxella osloensis",
                  "Staphylococcus epidermidis","Corynebacterium tuberculostearicum",
                  "Streptococcus oralis","Homo sapiens")
# PLOT
testing_confidenceSpec%>%
  mutate(truePos=if_else(name %in% true_positives,
                         name ,"False Positive"))%>%
  ggplot(aes(x=X8,y=fct_reorder(.x=X8,.f = truePos,.fun = median),fill= truePos!="False Positive"))+
  ggridges::geom_density_ridges(alpha=0.5)+
  scale_x_continuous(breaks = scales::pretty_breaks(n=10),limits = c(0,1))+
  geom_vline(xintercept = 0.1,lty=2)+
  labs(y="",x="Assignment Confidence")+
  guides(fill="none")+
  annotate(geom="text",x=0.1,y=13.5,label="Threshold=0.1",hjust=-.05,family="Roboto Condensed")
ggsave("plots/current_plots/confidences_benchmarking.pdf",device = cairo_pdf,height=5,width=7)
ggsave("plots/current_plots/confidences_benchmarking.png",device = png,height=5,width=7)

#### Supplementary Figure 9 - Contamination fractions in 16S data using SCRuB
# genus-level composition
negatives_composition<-
  s16_kraken%>%
  select(tax_id,starts_with("NTC"),starts_with("Water"))%>%
  filter(tax_id %in% (unified_taxonomy%>%
                        filter(type=="genus",str_detect(taxonomy,"Bacteria"))%>%
                        pull(tax_id)))%>%
  pivot_longer(-1, names_to = "Barcode",values_to = "Reads")%>%
  mutate(Relabund=Reads/sum(Reads),.by=Barcode)%>%
  left_join(s16_metadata,by=c("Barcode"="SampleID"))%>%
  left_join(unified_taxonomy)%>%
  graph_phyla_composition_barplot2(level = "genus", ordering_taxa = "Other", otus = 20) +
  facet_grid2(~ `Source Material` ,
              scales = "free",
              independent = "x",
              render_empty = F)+
  xlab(NULL)

# merge per-sample contamination fraction values
contamination_percents<-
  lapply(scrub_results_full,
         function(x){
           lapply(x$inner_iterations,function(x){1-x$p})%>%do.call(what="rbind")
         })%>%
  do.call(what="cbind")%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column("SampleID")
contamination_percents_simple<-
  lapply(scrub_results_full,
         function(x){1-x$p})%>%
  `names<-`(NULL)%>%
  unlist()%>%
  as.data.frame()%>%
  rownames_to_column("SampleID")%>%
  rename(Proportion=".")

# plot contamination percents
contamination_percents_plot<-
  contamination_percents%>%
  pivot_longer(-1,names_to = "Control",values_to = "Contamination Fraction")%>%
  left_join(s16_metadata)%>%
  ggplot(aes(x=PCR_Plate,y=`Contamination Fraction`,color=Control))+
  geom_boxplot(aes(fill=NULL),outlier.shape = NA)+
  geom_quasirandom(alpha=0.2,dodge.width = 0.8)

plot_grid(negatives_composition,
          contamination_percents_plot,
          ncol=1,labels = c("a","b"))
ggsave("plots/current_plots/contamination_fractions.png",device = png,width=7,height=7)
ggsave("plots/current_plots/contamination_fractions.pdf",device = cairo_pdf,width=7,height=7)