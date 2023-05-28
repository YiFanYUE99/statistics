install.packages("vegan")
library(ggplot2)
library(vegan)
#abundance: 
getwd()
#列名为样本，行名为代谢物丰度
metabolites<- read.csv("代谢物-数据文件.csv", row.names = 1, check.names = FALSE)
metabolite<-t(metabolites)
group<- data.frame(rownames(metabolite))# copy rownames

group$site<-c(rep('C',60),rep('P',60))#为样本分组

#使用PERMANOVA分析：Bray-Curtis距离，置换检验999次
adonis_result<-adonis2(metabolite~site,group,distance='bray',permutations=999)
adonis_result
#P=0.001

#输出
setwd("permanova")
getwd()
output<- data.frame(adonis_result,check.names = FALSE,stringsAsFactors = FALSE)
write.table(output,file = 'PERMANOVA_RESULT.txt',row.names = TRUE,sep = "\t",quote = FALSE)

#提取R2与P值
metabolite_adonis_R2<-paste("R2 :", round(adonis_result$R2[1],2))
metabolite_adonis_P<-paste("P-value :", round(adonis_result$`Pr(>F)`[1],2))

#做PCoA图
#PCoA作图
bray_dis <- vegdist(metabolite, method = 'bray')    
pcoa <- cmdscale(bray_dis, k = 2, eig = TRUE)
pcoa_eig <- pcoa$eig
pcoa_exp <- pcoa$eig/sum(pcoa$eig)
site <- scores(pcoa)

pcoa1 <- paste('PC1（', round(100*pcoa_exp[1], 2), '%）')
pcoa2 <- paste('PC2（', round(100*pcoa_exp[2], 2), '%）')

site <- data.frame(pcoa$point)[1:2]
site$name <- rownames(site)
site$treat <- c(rep('C', 60), rep('P', 60))

p <- ggplot(data = site, aes(X1, X2)) +
  geom_point(aes(color = treat)) +
  stat_ellipse(aes(fill = treat), geom = 'polygon', level = 0.95, alpha = 0.5, show.legend = FALSE) +   #添加置信椭圆，注意不是聚类
  scale_color_manual(values = c( 'green3', 'pink3')) +#change
  scale_fill_manual(values = c( 'green3', 'pink3')) +#change
  theme(panel.grid.major = element_line(),panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5), legend.position = "right") +
  geom_vline(xintercept = 0, color = 'gray', size = .1) +
  geom_hline(yintercept = 0, color = 'gray', size = .1) +
  labs(x = pcoa1, y = pcoa2, title = 'PCoA(PATIENTS_metabolite)')+
  annotate('text', label = metabolite_adonis_R2, x = -0.1, y = 0.23, size = 5, colour = 'black')+#change
  annotate('text', label = metabolite_adonis_P, x = -0.1, y = 0.2, size = 5, colour = 'black')#change
p

#两两比对
library(pairwiseAdonis)
#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
metabolite.pairwise.adonis <- pairwise.adonis(x=metabolite, factors=group$site, sim.function = "vegdist",
                                       sim.method = "bray",
                                       p.adjust.m = "BH",
                                       reduce = NULL,
                                       perm = 999)

metabolite.pairwise.adonis
metabolite.pairwise.adonis <- data.frame(metabolite.pairwise.adonis, stringsAsFactors = FALSE)
write.table(metabolite.pairwise.adonis, 'metabolite.pairwise.adonis.txt', row.names = TRUE, sep = '\t', quote = FALSE, na = '')

#在图中输出相关数据
library(ggpubr)
library(patchwork)
#输出表格ggtexttable
table <- ggtexttable(metabolite.pairwise.adonis[,c("pairs","R2","p.value","p.adjusted")], rows = NULL, 
                     theme = ttheme("blank")) %>% 
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 1)  %>% 
  tab_add_hline(at.row = nrow(metabolite.pairwise.adonis)+1, row.side = "bottom", linewidth = 1)  
p2 = p + table 
p2
#更改布局
patient_metabolite<-p2 + plot_layout(design=c(area(1,1), area(2,1)))
patient_metabolite
#输出图片到当前工作目录下
title_name = "patient_metabolite"
png(filename=paste0(title_name ,".png"),#好的批量命名方式paste0
    width=700,height=500,units="px",
    pointsize = 20)#标签大小
ggsave("patient_metabolite.png")
