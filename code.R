
library(dplyr)
install.packages("clusterProfiler")
library(clusterProfiler)

SNPlist <- read.xlsx("SNPlist.xlsx",sheet = 1)
go <- clusterProfiler::read.gmt("c5.go.v2022.1.Hs.symbols.gmt")
CP <- clusterProfiler::read.gmt("c2.all.v2022.1.Hs.symbols.gmt.txt")

gopathway <- left_join(SNPlist,go,by="gene")
gopathway$term <- as.character(gopathway$term)
cppathway <- left_join(SNPlist,CP,by="gene")
cppathway <- cppathway[grep("^KEGG|^PID|^RWACTOME|^WIKIPATHWAYS|^BIOCARTA",cppathway$term),] 
cppathway$term <- as.character(cppathway$term)
x1 =   cppathway %>% group_split(term)  
score <- c()
  for (i in 1:length(x1)) {
    write.table(x1[[i]][,c(1,4,6)],file = paste0("CPSCORE/score",i),
                quote = F,row.names = F,col.names = F)
    score[i] <- paste0("score",i," ",x1[[i]][1,8]," ",x1[[i]][1,9])
    write.table(score,quote = F,row.names = F,file="CPSCORE/name_score.txt")
    write.table(unique(cppathway$SNP),quote = F,row.names = F,file="CPSCORE/SNPlist.txt")
  }

x1 =   gopathway %>% group_split(term)  
score <- c()
for (i in 1:length(x1)) {
  write.table(x1[[i]][,c(1,4,6)],file = paste0("GOSCORE/score",i),
              quote = F,row.names = F,col.names = F)
  score[i] <- paste0("score",i," ",x1[[i]][1,8]," ",x1[[i]][1,9])
  write.table(score,quote = F,row.names = F,file="GOSCORE/name_score.txt")
  write.table(unique(cppathway$SNP),quote = F,row.names = F,file="CPSCORE/SNPlist.txt")
}
