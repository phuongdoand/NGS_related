fi = list.files(pattern = ".csv")
fi = fi[-grep("\\w{9}",fi)]
tar = read.csv("miRTarBase_MTI.csv")
tar_h = tar[which(tar[,3]=="Homo sapiens"),]

for(n in 1:6){
  DEmir = read.csv(fi[n])
  DEmir_u = DEmir[which(DEmir[,7]<0.05&DEmir[,3]>1),]
  DEmir_d = DEmir[which(DEmir[,7]<0.05&DEmir[,3]<(-1)),]
  tar_fu = tar_h[which(tar_h[,2]%in%DEmir_u[,1]),c(2,4)]
  tar_fd = tar_h[which(tar_h[,2]%in%DEmir_d[,1]),c(2,4)]
  write.csv(tar_fu,paste0(gsub(".csv","",fi[n]),"_",strsplit(fi[n],"")[[1]][1],"up_miRT.csv"),row.names = F)
  write.csv(tar_fd,paste0(gsub(".csv","",fi[n]),"_",strsplit(fi[n],"")[[1]][4],"up_miRT.csv"),row.names = F)
}

