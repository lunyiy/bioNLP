up = read.table(file = "G:\\Program Files\\data\\R\\GO.txt", sep = '\t', header = T, quote = "")
up_rt = up[up$PValue < 0.05,]

install.packages("tidyr")
library(tidyr)
up_rt = separate(up_rt, Term, sep = "~", into = c("ID", "Term"))

bp_df = up_rt[up_rt$Category == 'GOTERM_BP_DIRECT',]
bp_df = bp_df[order(bp_df$Count, decreasing = T),]
bp = bp_df[1:5,]

cc_df = up_rt[up_rt$Category == 'GOTERM_CC_DIRECT',]
cc_df = cc_df[order(cc_df$Count, decreasing = T),]
cc = cc_df[1:5,]

mf_df = up_rt[up_rt$Category == 'GOTERM_MF_DIRECT',]
mf_df = mf_df[order(mf_df$Count, decreasing = T),]
mf = mf_df[1:5,]

allGo = rbind(bp,cc,mf)

#install.packages("stringr")
library(stringr)
table(allGo$Category)
allGo$Category = substr(allGo$Category,8,9)

#install.packages("ggpubr")
library(ggplot2)
library(ggpubr)
colnames(allGo)
p = ggbarplot(data = allGo,x = "ID",y = 'Count',
              fill = "Category",
              palette = c("cadetblue3","mediumslateblue","mediumorchid3"),
              sort.by.groups = T,xlab = '',ylab = "Target genes")
ggpar(p,x.text.angle = 90)
ggsave(plot = p,'barplot.pdf',width = 10,height = 4)

