# This code was used to visualize the long-read sequencing data (figure 2A, supp. fig. 2, supp. fig. 5A). 

library(DNAcopy)
library(ggplot2)
library(reshape2)
library(scales)
library(grid)
library(gridExtra)
library(cowplot)

txtFontSize=16
axisFontSize=16
axisTtlFontSize=16
lgdTtlFontSize=16
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), 
                   legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), 
                   strip.background = element_blank(), axis.line=element_line(linewidth=0.7, color="black"), 
                   axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), 
                   axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), 
                   legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), 
                   text=element_text(size=txtFontSize), strip.text.x=element_text(size=axisTtlFontSize))

fmt_decimals = function(decimals=2) {
    function(x) format(x, nsmall=decimals, scientific=F)
    }

gg_color_hue = function(n) {
    hues=seq(15,375,length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
}

### figure 2A, supplementary figure 2 (all chromosome coverage by sample)
# uses hg38-aligned coverage

# replace with other samples for supplementary figure 2
sample <- "B123"

# load coverage file
cov_path <- paste0(sample, "tumor.cov")

plot_list <- list()
widths <- c()

minCN = 0
maxCN = 28

for (chrlabel in unique(cov$chr)) {
    
    # pull chromosome coverage
    reschr = cov[cov$chr==chrlabel,]$cn
    mid = (cov[cov$chr==chrlabel,2] + cov[cov$chr==chrlabel,3])/2
    
    df = data.frame(pos=mid, signal=reschr)
    
    # calculate width for this chromosome
    chr_width <- max(mid) - min(mid)
    widths <- c(widths, chr_width)
    
    # create facet plot for this chromosome
    p1 = ggplot(data=df, aes(x=pos, y=signal)) + geom_point(size=0.5, pch=21, fill="black", color="black")
    p1 = p1 + xlab(chrlabel) + ylab(NULL)
    p1 = p1 + scienceTheme
    p1 = p1 + scale_x_continuous(labels=function(x) paste0(x/1e6, "Mb"), limits=c(min(df$pos, na.rm=TRUE), max(df$pos, na.rm=TRUE)))
    p1 = p1 + scale_y_continuous(labels=NULL, breaks = seq(minCN, maxCN, by=4), limits=c(minCN, maxCN))
    p1 = p1 + theme(legend.position="bottom", axis.text.x = element_text(angle = 45, hjust = 1))
    
    plot_list[[chrlabel]] <- p1
}
                                 
rel_widths <- widths / sum(widths)
final_plot <- plot_grid(plotlist = plot_list, ncol = length(plot_list), align = "h", rel_widths = rel_widths)
                                 
options(repr.plot.width = 60)
options(repr.plot.height = 4)

print(final_plot)
                                 
### supplementary figure 5A (coverage for chromosome Y)  
# uses t2t-aligned coverage

# replace with other samples                              
sample = "B123"
cov_path <- paste0(sample, "tumor.cov")
idname = sample
cov = read.table(cov_path, header=T)
cov = cov[,c(1,2,3,6)]
colnames(cov) = c("chr", "start", "end", "cn")

plot_list <- list()
widths <- c()

chrlabel <- "chrY"
    
reschr = cov[cov$chr == chrlabel, ]$cn
mid = (cov[cov$chr == chrlabel, 2] + cov[cov$chr == chrlabel, 3]) / 2

df = data.frame(pos = mid, signal = reschr)

# compute average CN
avg_cn <- mean(reschr, na.rm = TRUE)

minCN = 0
maxCN = 4

chrY_plot = ggplot(data = df, aes(x = pos, y = signal)) +
  geom_point(size = 0.5, pch = 21, fill = "black", color = "black") +
  geom_hline(yintercept = avg_cn, linetype = "dashed", color = "red") +
  xlab(chrlabel) +
  ylab("Copy Number") +
  scale_x_continuous(
    labels = function(x) paste0(x / 1e6, "Mb"),
    limits = c(min(df$pos, na.rm = TRUE), max(df$pos, na.rm = TRUE))
  ) +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.5),
    breaks = seq(minCN, maxCN, by = 1),
    limits = c(minCN, maxCN)
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
                                 
options(repr.plot.width = 3)
options(repr.plot.height = 3)
print(chrY_plot)