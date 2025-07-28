# This code was used to plot L1 count versus SV count (figure 3F).

library(ggplot2)
library(dplyr)

txtFontSize=16
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16

scienceTheme=theme(panel.grid.major=element_blank(), 
                   panel.grid.minor=element_blank(), 
                   legend.key=element_blank(), 
                   legend.background=element_blank(), 
                   panel.background = element_blank(), 
                   panel.border=element_blank(),
                   strip.background = element_blank(), 
                   axis.line=element_line(linewidth=0.7, color="black"), 
                   axis.text.x=element_text(size=axisFontSize), 
                   axis.text.y=element_text(size=axisFontSize), 
                   axis.title.x=element_text(size=axisTtlFontSize), 
                   axis.title.y=element_text(size=axisTtlFontSize,angle = 90), 
                   legend.title=element_text(size=lgdTtlFontSize, 
                                             face="bold"),
                   legend.text=element_text(size=lgdFontSize),
                   text=element_text(size=txtFontSize), 
                   strip.text.x=element_text(size=axisTtlFontSize))
theme_set(scienceTheme)

df <- read.csv("l1_amp_sv_table.csv")

df <- df %>%
  mutate(
    amplicon_label = case_when(
      amplicon_high == 'True'  ~ "high (≥ 20)",
      amplicon_high == 'False' ~ "low (≤ 5)"
    )
  )

spearman_cor <- cor.test(df$sv_count, df$l1_count_log, method = "spearman")

p <- ggplot(df, aes(x = sv_count, y = l1_count_log, color = amplicon_label)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "gray80", linewidth = 1.2) +
  scale_color_manual(values = c("high (≥ 20)" = "purple", "low (≤ 5)" = "orange")) +
  labs(
    x = "Non-insertion SV count",
    y = "log(L1 count + 1)",
    title = "SV count versus L1 count",
    color = "Amplicon count"
  ) +
  annotate("text",
           x = min(df$sv_count, na.rm = TRUE) + 1,
           y = max(df$l1_count_log, na.rm = TRUE),
           label = sprintf("Spearman r = %.2f\np = %.3f", spearman_cor$estimate, spearman_cor$p.value),
           hjust = 0, vjust = 1, size = 5,
           fontface = "plain",
           color = "black",
           label.size = 0.2,
           fill = "white") +
  theme(
    legend.position = c(0.75, 0.2),
    legend.box.background = element_rect(color = "black", linewidth = 0.5),
    legend.background = element_blank()
  )

options(repr.plot.height = 4, repr.plot.width = 6)

print(p)