library(ggplot2)
library(scales)

txtFontSize=16
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), strip.text.x=element_text(size=axisTtlFontSize))


args=commandArgs(trailingOnly=TRUE)
x=read.table(args[1], header=T, sep="\t")
x$offset = x$pos
x[2:nrow(x),]$offset = x[2:nrow(x),]$pos - x[1:(nrow(x)-1),]$pos
x[1,]$offset = median(x$offset)

chrlabel = unique(x$chr)
print(args[2])
print(chrlabel)
cutoff = median(x$offset)/15
rleOff = rle(x$offset < cutoff)
x$rle = rep(rleOff$lengths, rleOff$lengths)
x[!rep(rleOff$values, rleOff$lengths),]$rle = 0

start = 0
end = 0
df = data.frame()
for(i in 1:nrow(x)) {
     ## At least 12 mutations
      if (x[i,]$rle > 12) {
          if (start == 0) { start = x[i,]$pos; }
	  end = x[i,]$pos;
      } else {
          if (start != 0) { df = rbind(df, c(chrlabel, start, end, args[2])); }
	  start = 0
	  end = 0
      }
}
if (start != 0) { df = rbind(df, c(chrlabel, start, end, args[2])); }
if (nrow(df)) {
  colnames(df) = c("chr", "start", "end", "sample")
  df$start = as.numeric(as.character(df$start))
  df$end = as.numeric(as.character(df$end))
  write.table(df, file=paste0("kataegis.", args[2], ".", chrlabel, ".tsv"), sep="\t", quote=FALSE, row.names=FALSE)
}

png(paste0("plot.", args[2], ".", chrlabel, ".png"), width=1200, height=600)
p = ggplot(data=x, aes(x=pos, y=offset))
p = p + geom_point(aes(color=mutation))
p = p + scale_y_log10(labels=comma)
p = p + scale_x_continuous(labels=comma)
p = p + ylab("Offset to prev. mutation (in log-scale)")
p = p + xlab(chrlabel)
p = p + scienceTheme
p = p + ggtitle(args[2])
if (nrow(df)) {
 p = p + geom_vline(data=df, aes(xintercept=start), linetype="dashed", color="lightgrey")
 p = p + geom_vline(data=df, aes(xintercept=end), linetype="dashed", color="lightgrey")
}
p
dev.off()



