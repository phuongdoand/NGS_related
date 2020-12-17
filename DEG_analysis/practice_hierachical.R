library(R.basic)
library(gplots)
library(ggplot2)
library(factoextra)
library(ggfortify)
library(RColorBrewer)

data = read.table('GSE26566_series_matrix.txt', skip = 60, header = T, nrows = 22184)
data = data[1:4000,]
row.names(data) = data[,1]
data = data[,2:(ncol(data)-6)]

# Just followed WeiNing, need to ask!!
data[data<0] = NaN

for(i in c(1:ncol(data))){
  a = which(is.na(data[,i]))
  data[a,i] = min(data[,i],na.rm = T)
}

t_z_score = t(apply(log2(data),1,zscore))
bk = c(seq(-2,-0.0011,by = 0.1),-0.001,0.001,seq(0.0011,2,by=0.1))
mycols = colorRampPalette(colors = c("blue","white","red"))(length(bk)-1)
my.hclust = function(d){
  hclust(d,method = 'ward.D')}

aa = heatmap.2(as.matrix(t_z_score), key = FALSE,
              trace = "none", breaks = bk, margins = c(1,1),
              col = mycols, scale = "row", Rowv = T,
              Colv = T, hclustfun = my.hclust, main = "GSE26566",
              cexRow = 1.5)

k = 6

cut_tree<-hcut(t(t_z_score),k = k,hc_func = 'hclust',hc_method = 'ward.D', hc_metric = 'euclidean')
cluster <- as.matrix(cut_tree$cluster)
cluster<-apply(cluster, 1, toString)
mcpscore_group <- cbind(t(t_z_score),cluster)
autoplot(prcomp(t(t_z_score)), main = paste0("Gse26566 euclidean & ward.D ", k," clusters"),
         data = mcpscore_group,colour = "cluster",size = 2)

my.dist<-function(x){ 
  as.dist(1-cor(t(x), method= dfun))}

d<-dist(t(t_z_score))
h<-my.hclust(d)
plot(h)
col1 <- brewer.pal(12, "Set3")

aa<-heatmap.2(as.matrix(t_z_score),key = FALSE,trace = "none",
              breaks = bk,col = mycols,scale = "row",
              Rowv = T,
              Colv = T,
              margins = c(0.2,1),
              hclustfun = my.hclust,
              main = 'GSE26566',
              # cexRow = 3.5,
              ColSideColors=col1[as.numeric(cluster)]
)