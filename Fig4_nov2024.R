source("RNA_survival.R")
source("WGS_survival.R")
source("16S_survival.R")
source("Survival_adenos_only.R")

# FOREST PLOT
survival_top<-
  plot_grid(genera_survival_plot+
              guides(color="none")+
              theme(axis.text.x = element_text(size=8)),
            genera_16S_survival_plot+
              guides(color="none")+
              theme(axis.text.x = element_text(size=8)),
            meta_wgs_survival_plot+
              guides(color="none")+
              theme(axis.text.x = element_text(size=8)),
            ncol=3,
            labels=c("a","b","c"),
            align="hv")
# COLOR LABELS
survival_bottom<-
  plot_grid(NULL,
            get_legend(genera_16S_survival_plot+theme(legend.direction = "horizontal")),
            NULL,
            get_legend(meta_wgs_survival_plot+theme(legend.direction = "horizontal")),
            nrow=1)
plot_grid(survival_top,
          survival_bottom,
          ncol=1,
          rel_heights = c(1,.075))
ggsave("plots/current_plots/survival_genera_v2.pdf",device=cairo_pdf,height=9,width=10,bg = "#ffffff")
ggsave("plots/current_plots/survival_genera_v2.png",height=9,width=10,bg = "#ffffff")