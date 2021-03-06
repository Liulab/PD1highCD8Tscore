getSigScore <-
function(tpm,signature_genes,background=background){
rank_list<-rankGenes(tpm[intersect(row.names(tpm),background),], tiesMethod = "min")
if (length(signature_genes)==1){
score<-simpleScore(rank_list,signature_genes[[1]],subSamples = NULL, centerScore =F,dispersionFun = mad)
}else{
score<-simpleScore(rank_list,upSet=signature_genes[[1]],downSet=signature_genes[[2]],subSamples = NULL, centerScore =F,dispersionFun = mad)
}
gscore<-score[,1]
names(gscore)<-row.names(score)
return(gscore)
}
