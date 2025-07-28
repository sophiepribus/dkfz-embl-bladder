# This code was used to plot Y signature score versus average chrY CN (supplementary fig. 5b).

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

df <- read.csv("072825_chrYCN_score.csv")

spearman_cor <- cor.test(df$avg_chrY_CN, df$Y_score, method = "spearman")

p <- ggplot(df, aes(x = avg_chrY_CN, y = Y_score)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "gray80", linewidth = 1.2) +
  labs(
    x = "Average chrY copy number",
    y = "Y signature score",
    title = "Y signature score vs average chrY CN"
  ) +
  annotate("label",
           x = min(df$avg_chrY_CN, na.rm = TRUE),
           y = max(df$Y_score + 10000, na.rm = TRUE),
           label = sprintf("Spearman r = %.2f\np = %.3f", spearman_cor$estimate, spearman_cor$p.value),
           hjust = 0, vjust = 1, size = 5,
           fontface = "plain",
           color = "black",
           label.size = 0.2,
           fill = "white")

options(repr.plot.height = 4, repr.plot.width = 6)

print(p)