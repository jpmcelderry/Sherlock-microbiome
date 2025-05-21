set.seed(123)
source("diversity_calculations.R")
source("RNA_diversity.R")
source("WGS_diversity.R")
source("16S_diversity.R")
source("RNA_differential_abundance.R")
source("16S_differential_abundance.R")

######### FIGURE 3b and 3d
diversity_data_RNA <-
  RNAseq_diversity %>%
  left_join(rna_annots) %>%
  filter(Type %in% c("Normal", "Tumor")) %>%
  pivot_longer(`500`:`1500`, names_to = "depth", values_to = "Shannon diversity") %>%
  filter(depth == "500") %>%
  drop_na(`Shannon diversity`) %>%
  filter(nlevels(as.factor(Type)) > 1,.by=c(Sherlock_PID, depth)) %>%
  # if there are multiple samples per subject+TN status, average
  summarise(`Shannon diversity`=mean(`Shannon diversity`),
            .by = c(Sherlock_PID,Type))
# N label for the plot
nSubject_aDiversity_rna<-
  nlevels(as.factor(diversity_data_RNA$Sherlock_PID))
# make violin plot
diversity_plot_RNA<-
  diversity_data_RNA%>%
  mutate(experiment = "RNA-seq")%>%
  ggplot(aes(x = Type, y = `Shannon diversity`)) +
  geom_violin() +
  geom_boxplot(width = 0.2, fill = NA, outlier.alpha = 0.4) +
  stat_compare_means() +
  labs(x = paste0("500 read cutoff, n=",nSubject_aDiversity_rna," pairs"))+
  ggh4x::facet_wrap2(~experiment,scales = "free_x",strip = strip_themed(background_x = element_rect(fill=RNA_color)))
diversity_plot_RNA

diversity_data_16S <-
  s16_diversity %>%
  pivot_longer(`100`:`750`,names_to = "depth",values_to = "Shannon diversity")%>%
  left_join(s16_clinical_annos)%>%
  filter(depth == "250", `Tumor-NormalStatus` %in% c("Tumor", "Normal")) %>%
  drop_na(`Shannon diversity`)%>%
  filter(nlevels(as.factor(`Tumor-NormalStatus`)) > 1,
         .by=AdditionalAttributes) %>%
  # if there are multiple samples per subject+TN status, average
  summarise(`Shannon diversity`=mean(`Shannon diversity`),
            .by = c(AdditionalAttributes,`Tumor-NormalStatus`))%>%
  rename(Type = "Tumor-NormalStatus")
# N label for the plot
nSubject_aDiversity_16S<-nlevels(as.factor(diversity_data_16S$AdditionalAttributes))
# make violin plot
diversity_plot_16S<-
  diversity_data_16S%>%
  mutate(experiment = "16S")%>%
  ggplot(aes(x = Type, y = `Shannon diversity`)) +
  geom_violin() +
  # geom_quasirandom(alpha=0.5)+
  geom_boxplot(width = 0.2, fill = NA, outlier.alpha = 0.4) +
  stat_compare_means() +
  labs(x = paste0("250 read cutoff, n=",nSubject_aDiversity_16S," pairs"))+
  ggh4x::facet_wrap2(~experiment,scales = "free_x",
                     strip = strip_themed(background_x = element_rect(fill=s16_color)))
diversity_plot_16S

multiOme_forest<-
  wgs_richness_diversity_results%>%
  rbind(s16_richness_diversity_results%>%select(-statistic))%>%
  rbind(RNA_richness_diversity_results%>%select(-statistic))%>%
  mutate_at("p.value",round,digits=2)%>%
  arrange(term)%>%
  forestplot(name = term,
             estimate = estimate,
             se=std.error)+
  labs(x="Correlation coefficient")+
  geom_vline(xintercept = 0,lty=3,alpha=0.5)+
  geom_text(mapping = aes(x=estimate,y=term,label=paste0("p=",round(p.value,digits = 3))),
            nudge_y = 0.25,family="Roboto Condensed",size=2.5)+
  theme_bw(base_family = "Roboto Condensed")+
  theme(text = element_text(size = 10))+
  ggh4x::facet_wrap2(~var,scales = "free_x",nrow = 1,
                     strip = strip_themed(
                       background_x = elem_list_rect(
                         fill=c(s16_color,s16_color,
                                RNA_color,RNA_color,
                                WGS_color,WGS_color))))
multiOme_forest$layers[[2]]<-NULL
multiOme_forest

f3_top <- plot_grid(differential_abundance_16S + 
                      theme(legend.position = "none"),
                    diversity_plot_16S,
                    differential_abundance_RNA,
                    diversity_plot_RNA,
                    labels = c("a", "b", "c", "d"),
                    nrow = 1, axis = "b", align = "h",
                    rel_widths = c(1, 1, 1, 1))
plot_grid(f3_top, multiOme_forest, 
          nrow = 2,rel_heights = c(0.45,0.55),labels=c("","e"))

ggsave("plots/current_plots/nonfindings_figure_v2.pdf",device=cairo_pdf,height = 10,width = 11)
ggsave("plots/current_plots/nonfindings_figure_v2.png",height = 10,width = 11)