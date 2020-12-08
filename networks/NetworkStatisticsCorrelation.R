setwd("/Users/mguardado1/Desktop/UCSF/Coursework/BMI206/FinalProject/networks/")
library(reshape2)
library(ggplot2)
library(caret)
library(klaR)


K_df=read.csv("OrganismGraphInformation.csv")

K_orgnames=unique(K_df$organism)

for (i in K_orgnames){
  TempDf=subset(K_df, organism==i & k_type=="total_k")
  print(i)
  print(summary(lm(log_p~log_k,data=TempDf)))
  
}

summary(lm(log_p~log_k,data=K_df))

NetworkStatistics=read.csv("MetabolomicsNetworkStatistics.csv")
NetworkStatistics$X<-NULL
NetworkStatistics$k_prob<-NULL

library(GGally)
ggpairs(NetworkStatistics[2:9]) +theme_classic()

NetworkCorr=round(cor(as.matrix(NetworkStatistics[2:9])),2)
MeltedNetworkCorr=melt(NetworkCorr)


ggplot(data = MeltedNetworkCorr, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

Network_upper_tri=get_upper_tri(NetworkCorr)

melted_cormat <- melt(Network_upper_tri, na.rm = TRUE)
# Heatmap
library(ggplot2)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

dglm.fit=lm(CascadeNum~ k_in+k_out+BetweennessCentrality+ClusteringCoeff+Bridge+BridgingCentrality,data=NetworkStatistics)


NetworkStatSubset=subset(NetworkStatistics, CascadeNum != 0)
NetworkStatSubsetCorr=round(cor(as.matrix(NetworkStatSubset[2:9])),2)

NetworkStatSubsetUpperTriangle=get_upper_tri(NetworkStatSubsetCorr)

melted_cormat <- melt(NetworkStatSubsetUpperTriangle, na.rm = TRUE)
# Heatmap
library(ggplot2)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()


cascadeFit = lm(CascadeNum~ k_total+BetweennessCentrality+BridgingCentrality+Organism, data=NetworkStatistics)

print(summary(cascadeFit))
unique(NetworkStatistics$Organism)







library(reshape2)
library(ggplot2)
library(caret)
library(klaR)
NetworkStatistics=read.csv("MetabolomicsNetworkStatistics.csv")
NetworkStatistics$X<-NULL


NetworkStatisticsSub=NetworkStatistics[c('CascadeNum','k_total','BetweennessCentrality','BridgingCentrality')]
#This will not create a model to test and trian the model to test the accuracy of predictice cascade number
set.seed(263)
smp_size <- floor(0.80 * nrow(NetworkStatistics))

## set the seed to make your partition reproducible
train_ind <- sample(seq_len(nrow(NetworkStatisticsSub)), size = smp_size)

train <- NetworkStatisticsSub[train_ind, ]
test <- NetworkStatisticsSub[-train_ind, ]

fit.1 = lm(CascadeNum~ k_total+BetweennessCentrality+BridgingCentrality, data=train)
fit.2 = lm(CascadeNum~ k_total*BetweennessCentrality*BridgingCentrality, data=train)

Predfit.1 <- round(predict(fit.1, test),0)
Predfit.2 <- round(predict(fit.2, test),0)

actuals_Predfit.1 <- data.frame(cbind(actuals=test$CascadeNum, predicteds=Predfit.1))  # make actuals_predicteds dataframe.

length(which(Predfit.1==test$CascadeNum))/length(test$CascadeNum)

actuals_Predfit.2 <- data.frame(cbind(actuals=test$CascadeNum, predicteds=Predfit.2))  # make actuals_predicteds dataframe.

length(which(Predfit.2==test$CascadeNum))/length(test$CascadeNum)

#Fitting a GLM model for a poison distribution

glmfit.1 = glm(CascadeNum~ k_total+BetweennessCentrality+BridgingCentrality, 
           family = poisson(), data=train)

glmfit.2 = glm(CascadeNum~ k_total*BetweennessCentrality*BridgingCentrality, 
           family = poisson(link = "log"),data=train)

summary(glmfit.1)
summary(glmfit.2)

Predfit.1 <- round((predict(glmfit.1, test,type = "response")),0)

Predfit.2 <- round(predict(glmfit.1, test,type = "response"),0)

truth<-as.factor(as.numeric(test$CascadeNum))

actuals_Predfit.1 <- data.frame(cbind(actuals=truth, predicteds=Predfit.1))
# make actuals_predicteds dataframe.
length(which(Predfit.1==truth))/length(truth)

AIC(glmfit.1)

actuals_Predfit.2 <- data.frame(cbind(actuals=test$CascadeNum, predicteds=Predfit.2))  # make actuals_predicteds dataframe.
length(which(Predfit.2==truth))/length(truth)


glmfit.3 = glm.nb(CascadeNum ~ k_total+BetweennessCentrality+BridgingCentrality, data=train)



