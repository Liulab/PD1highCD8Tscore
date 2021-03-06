getMeanExpression <-
function(EXP,x){
x<-intersect(x,row.names(EXP))
house<-c('RPL38', 'UBA52', 'RPL4', 'RPS29', 'SLC25A3', 'CLTC', 'RPL37', 'PSMA1', 'RPL8', 'PPP2CA','TXNL1', 'MMADHC', 'PSMC1', 'RPL13A', 'MRFAP1')
house<- apply(EXP[house,],2,mean)
if (length(x)==1){
ME<-EXP[x,]
}else{
ME<-apply(EXP[x,],2,mean)}
ME<-ME/house
return(ME)
}
