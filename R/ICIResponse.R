ICIResponse <-
function(marker,Gide=Gide,Riaz=Riaz,Kim=Kim,Cho=Cho,Jung=Jung,background=background){
DF<-list()
if (class(marker)=='character'){
marker_Gide<-getMeanExpression(Gide$RNA,marker)
marker_Riaz<-getMeanExpression(Riaz$RNA,marker)
marker_Kim<-getMeanExpression(Kim$RNA,marker)
marker_Jung<-getMeanExpression(Jung$RNA,marker)
marker_Cho<-getMeanExpression(Cho$RNA,marker)}
else{
marker_Gide<-getSigScore(Gide$RNA,marker,background=background)
marker_Riaz<-getSigScore(Riaz$RNA,marker,background=background)
marker_Kim<-getSigScore(Kim$RNA,marker,background=background)
marker_Jung<-getSigScore(Jung$RNA,marker,background=background)
marker_Cho<-getSigScore(Cho$RNA,marker,background=background)}

Gide_pre_No<-marker_Gide[Gide$clinic$Alias=='PRE' & Gide$clinic$'response'=='No']
Gide_pre_Yes<-marker_Gide[Gide$clinic$Alias=='PRE' & Gide$clinic$'response'=='Yes']
p1<-wilcox.test(Gide_pre_No,Gide_pre_Yes)$p.value
auc1<-roc(c(rep('No',length(Gide_pre_No)),rep('Yes',length(Gide_pre_Yes))),c(Gide_pre_No,Gide_pre_Yes))

Gide_On_No<-marker_Gide[Gide$clinic$Alias=='EDT' & Gide$clinic$'response'=='No']
Gide_On_Yes<-marker_Gide[Gide$clinic$Alias=='EDT' & Gide$clinic$'response'=='Yes']
p2<-wilcox.test(Gide_On_No,Gide_On_Yes)$p.value
auc2<-roc(c(rep('No',length(Gide_On_No)),rep('Yes',length(Gide_On_Yes))),c(Gide_On_No,Gide_On_Yes))

Riaz_pre_No<-marker_Riaz[Riaz$clinic$time=='Pre' & Riaz$clinic$response=='No']
Riaz_pre_Yes<-marker_Riaz[Riaz$clinic$time=='Pre' & Riaz$clinic$response=='Yes']
p3<-wilcox.test(Riaz_pre_No,Riaz_pre_Yes)$p.value
auc3<-roc(c(rep('No',length(Riaz_pre_No)),rep('Yes',length(Riaz_pre_Yes))),c(Riaz_pre_No,Riaz_pre_Yes))

Riaz_On_No<-marker_Riaz[Riaz$clinic$time=='On' & Riaz$clinic$response=='No']
Riaz_On_Yes<-marker_Riaz[Riaz$clinic$time=='On' & Riaz$clinic$response=='Yes']
p4<-wilcox.test(Riaz_On_No,Riaz_On_Yes)$p.value
auc4<-roc(c(rep('No',length(Riaz_On_No)),rep('Yes',length(Riaz_On_Yes))),c(Riaz_On_No,Riaz_On_Yes))

Kim_No<-marker_Kim[Kim$clinic$response=='No']
Kim_Yes<-marker_Kim[Kim$clinic$response=='Yes']
p5<-wilcox.test(Kim_No,Kim_Yes)$p.value
auc5<-roc(c(rep('No',length(Kim_No)),rep('Yes',length(Kim_Yes))),c(Kim_No,Kim_Yes))

Jung_No<-marker_Jung[Jung$clinic$response=='No']
Jung_Yes<-marker_Jung[Jung$clinic$response=='Yes']
p7<-wilcox.test(Jung_No,Jung_Yes)$p.value
auc7<-roc(c(rep('No',length(Jung_No)),rep('Yes',length(Jung_Yes))),c(Jung_No,Jung_Yes))

Cho_No<-marker_Cho[Cho$clinic$Clinicalbenefit=='non-responder']
Cho_Yes<-marker_Cho[Cho$clinic$Clinicalbenefit=='responder']
p8<-wilcox.test(Cho_No,Cho_Yes)$p.value
auc8<-roc(c(rep('No',length(Cho_No)),rep('Yes',length(Cho_Yes))),c(Cho_No,Cho_Yes))

DF[[1]]<-data.frame(DataSet=c('Gide_Pretreatment','Gide_Ontreatment','Riaz_Pretreatment','Riaz_Ontreatment','Kim','Jung','Cho'),
   No_response=c(mean(Gide_pre_No),mean(Gide_On_No),mean(Riaz_pre_No),mean(Riaz_On_No),mean(Kim_No),mean(Jung_No),mean(Cho_No)),
   response=c(mean(Gide_pre_Yes),mean(Gide_On_Yes),mean(Riaz_pre_Yes),mean(Riaz_On_Yes),mean(Kim_Yes),mean(Jung_Yes),mean(Cho_Yes)),
   WilCox_p=c(p1,p2,p3,p4,p5,p7,p8),
   AUC=c(as.numeric(auc1$auc),as.numeric(auc2$auc),as.numeric(auc3$auc),as.numeric(auc4$auc),as.numeric(auc5$auc),as.numeric(auc7$auc),as.numeric(auc8$auc)))

##Survival
Gide_pre<-Gide$clinic[Gide$clinic$Alias=='PRE',]
fit1<-coxph(Surv(Gide_pre$'Overall.Survival..Days.',Gide_pre$status!='Alive')~scale(marker_Gide[Gide$clinic$Alias=='PRE']))
fit1.1<-coxph(Surv(Gide_pre$'Progression.Free.Survival..Days.',Gide_pre$Progressed=='Yes')~scale(marker_Gide[Gide$clinic$Alias=='PRE']))

Gide_On<-Gide$clinic[Gide$clinic$Alias=='EDT',]
fit2<-coxph(Surv(Gide_On$'Overall.Survival..Days.',Gide_On$status!='Alive')~scale(marker_Gide[Gide$clinic$Alias=='EDT']))
fit2.1<-coxph(Surv(Gide_On$'Progression.Free.Survival..Days.',Gide_On$Progressed=='Yes')~scale(marker_Gide[Gide$clinic$Alias=='EDT']))


fit4<-coxph(Surv(Jung$clinic$PFS,Jung$clinic$PD_Event.1_Censoring.0)~scale(marker_Jung))

DF[[2]]<-rbind(coef(summary(fit1)),coef(summary(fit1.1)),coef(summary(fit2)),coef(summary(fit2.1)),coef(summary(fit4)))
row.names(DF[[2]])<-c('Gide_Pre_OS','Gide_Pre_PFS','Gide_On_OS','Gide_On_PFS','Jung_PFS')

names(DF) <- c("AUC","cox")
return(DF)}
