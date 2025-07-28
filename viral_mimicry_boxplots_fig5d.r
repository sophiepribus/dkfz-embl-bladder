# This code was used for the viral mimicry pathway analysis (figure 5D). 

library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(stats)

txtFontSize=16
axisFontSize=16
axisTtlFontSize=18
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

# load gene expression data
all_genes_df <- read_tsv("060225_salmon.merged.gene_counts.tsv")

# define sample groups
l1_high <- c("B42tumor", "B123tumor", "B134tumor", "B4tumor", "B5tumor")
l1_low <- c("B24tumor", "B154tumor", "B178tumor", "B39tumor", "B60tumor", "B156tumor", "B74tumor")

plot_gene_boxplot <- function(df, gene) {
  df_long <- df %>%
    filter(gene_name == gene) %>%
    select(-gene_id) %>%
    pivot_longer(-gene_name, names_to = "sample", values_to = "expr") %>%
    mutate(
      expr = as.numeric(expr),
      log_expr = log10(expr + 1),
      L1_status = case_when(
        sample %in% l1_high ~ "L1-high",
        sample %in% l1_low ~ "L1-low",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(L1_status))

  # Mann-Whitney U test
  x <- df_long$log_expr[df_long$L1_status == "L1-high"]
  y <- df_long$log_expr[df_long$L1_status == "L1-low"]
  test <- wilcox.test(x, y, alternative = "greater", exact = FALSE)

  n1 <- length(x)
  n2 <- length(y)
  u1 <- as.numeric(test$statistic)
  u2 <- n1 * n2 - u1

  # Adjusted rank-biserial effect size (with sign based on U1 vs U2)
  if (u1 > u2) {
    r_rb <- (u1 - u2) / (n1 * n2)
  } else {
    r_rb <- (u2 - u1) / (n1 * n2)
    r_rb <- -r_rb
  }
    
  print(test)
  cat("Rank-biserial effect size r:", round(r_rb, 4), "\n")

  p_val <- test$p.value
  label_text <- paste0("p = ", signif(p_val, 3))

  # Plot
  p <- ggplot(df_long, aes(x = L1_status, y = log_expr)) +
    geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +
    geom_jitter(width = 0.2, size = 4, color = "black", alpha = 0.6) +
    labs(
      title = gene,
      y = "log(gene count + 1)",
      x = NULL
    ) +
    scienceTheme +
    annotate(
      "text",
      x = 1.5,
      y = max(df_long$log_expr, na.rm = TRUE) * 1.1,
      label = label_text,
      size = 6
    )

  options(repr.plot.width=2.5, repr.plot.height=4)
  print(p)
}

plot_gene_boxplot(all_genes_df, "APOBEC3B")
plot_gene_boxplot(all_genes_df, "RIGI")
plot_gene_boxplot(all_genes_df, "IFIH1")

# define cGAS-STING pathway genes
cgas_sting_genes <- c("CGAS", "CCL5", "CXCL10", "IRF3", "TBK1", "STING1", "STAT1")

cgas_sting_df <- all_genes_df %>%
  filter(gene_name %in% cgas_sting_genes) %>%
  select(-gene_id)

cgas_sting_mat <- cgas_sting_df %>%
  column_to_rownames("gene_name") %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

# apply log2(x + 1) transformation
cgas_sting_log <- log2(cgas_sting_mat + 1)

# center
cgas_sting_centered <- t(scale(t(cgas_sting_log), center = TRUE, scale = FALSE))

# compute cGAS-STING scores (sum of centered expression)
cgas_sting_centered <- data.frame(
  sample = colnames(cgas_sting_centered),
  cGAS_STING_score = colSums(cgas_sting_centered, na.rm = TRUE)
)

cgas_sting_scores <- cgas_sting_centered %>%
  mutate(L1_status = case_when(
    sample %in% l1_high ~ "L1-high",
    sample %in% l1_low ~ "L1-low",
    TRUE ~ "other"
  )) %>%
  filter(L1_status %in% c("L1-high", "L1-low"))

plot_cgas_sting_boxplot <- function(df) {
  df_filtered <- df %>%
    filter(L1_status %in% c("L1-high", "L1-low"))

  # Mann-Whitney U-test
  x <- df_filtered$cGAS_STING_score[df_filtered$L1_status == "L1-high"]
  y <- df_filtered$cGAS_STING_score[df_filtered$L1_status == "L1-low"]
  test <- wilcox.test(x, y, alternative = "greater", exact = FALSE)

  n1 <- length(x)
  n2 <- length(y)
  u1 <- as.numeric(test$statistic)
  u2 <- n1 * n2 - u1

  # Adjusted rank-biserial effect size (with sign based on U1 vs U2)
  if (u1 > u2) {
    r_rb <- (u1 - u2) / (n1 * n2)
  } else {
    r_rb <- (u2 - u1) / (n1 * n2)
    r_rb <- -r_rb
  }
    
  print(test)
  cat("Rank-biserial effect size r:", round(r_rb, 4), "\n")

  p_val <- test$p.value
  label_text <- paste0("p = ", signif(p_val, 3))

  # Plot
  p <- ggplot(df_filtered, aes(x = L1_status, y = cGAS_STING_score)) +
    geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +
    geom_jitter(width = 0.2, size = 4, color = "black", alpha = 0.6) +
    labs(
      title = "cGAS-STING",
      y = "cGAS-STING score",
      x = NULL
    ) +
    scienceTheme +
    annotate("text", x = 1.5, y = max(df_filtered$cGAS_STING_score, na.rm = TRUE) * 1.15, label = label_text, size = 6)

  options(repr.plot.width=2.5, repr.plot.height=4)
  print(p)
}

plot_cgas_sting_boxplot(cgas_sting_scores)