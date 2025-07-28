# This code was used to visualize figure 5a.

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

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

vafs_df <- read.csv("072825_vafs.csv")

vafs_df <- vafs_df %>%
  mutate(SV.type = factor(SV.type), 
         Clonal = factor(Clonal))

proportions <- vafs_df %>%
  group_by(SV.type, Clonal) %>%
  tally() %>%
  group_by(SV.type) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup() %>%
  spread(key = Clonal, value = proportion, fill = 0) %>%
  # Collapse proportions for each SV type
  group_by(SV.type) %>%
  summarise(
    n = sum(n),              # Total count for each SV type
    False = sum(`False` * n) / sum(n),  # Weighted False proportion
    True = sum(`True` * n) / sum(n)     # Weighted True proportion
  ) %>%
  ungroup()

proportions <- proportions %>%
  arrange(desc(True))

proportions <- proportions %>%
  mutate(label = paste(SV.type, "\n(n=", n, ")", sep = ""))

df_long <- proportions %>%
  pivot_longer(cols = c(False, True), names_to = "status", values_to = "proportion")

# to remove SVAs and Alus
df_long <- df_long %>%
  filter(!(SV.type %in% c("INS - SVA", "INS - Alu")))

label_map <- df_long %>%
  distinct(SV.type, label) %>%
  deframe()

ordered_SVtype <- df_long %>%
  filter(status == "True") %>%
  arrange(desc(proportion)) %>%
  pull(SV.type)

options(repr.plot.width=10, repr.plot.height=6)
ggplot(df_long, aes(x = SV.type, y = proportion, fill = status)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("True" = "navy", "False" = "orange")) + 
  scale_x_discrete(
    limits = ordered_SVtype,
    labels = label_map
  ) +
  labs(y = "Proportion", fill = "Clonal", title = "Proportions of Clonal/Subclonal by SV Type") +
  theme_minimal() + 
  scienceTheme