setwd("/Users/mguardado1/Desktop/UCSF/Coursework/BMI206/FinalProject/networks/")
library(reshape2)
library(ggplot2)
library(tidyverse)
library(caret)
library(klaR)
library(plyr)
library(ggfortify)
library(PCAtools)
library(e1071)

#Load in Organism Data
NetworkStatistics=read.csv("MetabolomicsNetworkStatisticsUpdated.csv")
K_orgnames=unique(NetworkStatistics$Organism)
NetworkStatistics$X<-NULL
NetworkStatistics=NetworkStatistics[c('Organism','CascadeNum','k_total','BetweennessCentrality','BridgingCentrality')]

#Load in data based off the essentially of each node
EssentialReactionseColi=read.csv('ecoliES.csv')
EssentialReactionskPneu=read.csv('kpneumoniaeES.csv')
EssentialReactionsbSub=read.csv('bsubtilisES.csv')
EssentialReactionsgMeta=read.csv('gmetallireducensES.csv')
EssentialReactionsesCere=read.csv('scerevisiaeES.csv')

#Since I node the nodes in EssentialReactions are the same as the NetworkStatistics I will simply just 
#combine the two data frames together. 
Sub=subset(NetworkStatistics,Organism=='kPneumoniae')
Sub$Nodes=factor(Sub$Nodes)
all.equal(Sub$Nodes,EssentialReactionskPneu$node)

NetworkStatistics$EssentialReaction=c(EssentialReactionseColi$essential,EssentialReactionsbSub$essential,
                                      EssentialReactionsesCere$essential,EssentialReactionsgMeta$essential,
                                      EssentialReactionskPneu$essential)
#Lol my code is
NetworkStatistics$EssentialReaction <- mapvalues(NetworkStatistics$EssentialReaction, 
                                                 from = c(1, 2), to = c(0, 1))

#Create Test and Training data based on sampling on the the 5 Org as the 
set.seed(263)
TestOrg=sample(K_orgnames,1)
train <- subset(NetworkStatistics, Organism!=TestOrg)
test <- subset(NetworkStatistics, Organism==TestOrg)

#Redoing prediction of Cascade number based on updated testing and training datasets  
glmfit.1 = glm(CascadeNum~k_total+EigenCentrality+ClusteringCoeff, 
               family = poisson(), data=train)


summary(glmfit.1)

Predfit.1 <- round((predict(glmfit.1, test,type = "response")),0)
truth<-as.factor(as.numeric(test$CascadeNum))
length(which(Predfit.1==truth))/length(truth)
table(Predfit.1,truth)



#Now I will present an additional model that will look at a multinomial distribuition. i
library(nnet)
glmfit.2=multinom(CascadeNum~ k_total+BetweennessCentrality+BridgingCentrality, data=train)
print(summary(glmfit.2))

Predfit.2 <- predict(glmfit.2, newdata = test,type="class")
truth<-as.factor(as.numeric(test$CascadeNum))
length(which(Predfit.2==truth))/length(truth)

#######################################################################################################################

NetworkStatisticsSub=NetworkStatistics[c(2,3,4,5,6,7,8)]
NetworkStatisticsSub$k_total=as.numeric(NetworkStatisticsSub$k_total)
NetworkStatisticsSub$CascadeNum=as.numeric(NetworkStatisticsSub$CascadeNum)

pca_res <- prcomp(as.matrix(NetworkStatisticsSub), scale = TRUE, center=TRUE)
df_out <- as.data.frame(pca_res$x)
df_out$group <- NetworkStatistics$EssentialReaction
df_out$group <- mapvalues(df_out$group,  from = c(0, 1), to = c("Yes", "No"))

#pca_plot <- pca(NetworkStatisticsSub)
#screeplot(pca_res)
#plotloadings(pca_plot)


p1<-ggplot(df_out,aes(x=PC1,y=PC2,color=group))

p1<-p1+geom_point() + ggtitle('PC1 vs PC2')

p1
ggsave("PC1PC2.png")

p2<-ggplot(df_out,aes(x=PC2,y=PC3,color=group))

p2<-p2+geom_point() + ggtitle('PC2 vs PC3')

p2
ggsave("PC2PC3.png")


p3<-ggplot(df_out,aes(x=PC3,y=PC4,color=group))
p3<-p3+geom_point() + ggtitle('PC3 vs PC4')

p3
ggsave("PC3PC4.png")


#####SUB_PCA_PLOT >:)

NetworkStatisticsSub=NetworkStatistics[c(3,7)]
NetworkStatisticsSub$CascadeNum=as.numeric(NetworkStatisticsSub$CascadeNum)


pca_res <- prcomp(as.matrix(NetworkStatisticsSub), scale = TRUE, center=TRUE)
df_out <- as.data.frame(pca_res$x)

NetworkStatisticsSub$Essential <- NetworkStatistics$EssentialReaction
NetworkStatisticsSub$Essential <- mapvalues(NetworkStatisticsSub$Essential,  from = c(0, 1), to = c("Yes", "No"))



p4<-ggplot(df_out,aes(x=PC1,y=PC2,color=group))
p4<-p4+geom_point() + ggtitle('PCA of BridgingCentrality+CascadeNum')
p4
X




#######################################################################################################################
#This code will now look at trying to model predict

#Create Test and Training data based on sampling on the the 5 Org as the 
set.seed(263)
TestOrg=sample(K_orgnames,1)
train <- subset(NetworkStatistics, Organism!=TestOrg)
test <- subset(NetworkStatistics, Organism==TestOrg)

#Using just logistic regression
logfit.1 <- glm(EssentialReaction ~  CascadeNum * k_total* EigenCentrality, data = train, family = "binomial")
summary(logfit.1)

truth=test$EssentialReaction
truth

probabilities <- logfit.1 %>% predict(test, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
predicted.classes

Predfit.3 <- predict(logfit.1, newdata = test,type="response")

mean(predicted.classes == test$EssentialReaction)
table(predicted.classes,truth)


#Trying to Utilize Feature Selection 
library(glmnet)
x=as.matrix(NetworkStatisticsSub)
y=NetworkStatistics$EssentialReaction

cv.out <- cv.glmnet(x,y,alpha=1,family="binomial",type.measure = "mse" )
#plot result
plot(cv.out)

#min value of lambda
lambda_min <- cv.out$lambda.min
#best value of lambda
lambda_1se <- cv.out$lambda.1se

coef(cv.out,s=lambda_1se)


#Lasso regression found BetweenessCentrality and CascadeNum do the best job in predicitng Essential reactions 
logfit.1 <- glm(EssentialReaction ~  BetweennessCentrality+CascadeNum, data = train, family = "binomial")
summary(logfit.1)

truth=test$EssentialReaction
truth

probabilities <- logfit.1 %>% predict(test, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
predicted.classes

mean(predicted.classes == test$EssentialReaction)

##########
logfit.2 <- glm(EssentialReaction ~  BetweennessCentrality*CascadeNum, data = train, family = "binomial")
summary(logfit.2)

truth=test$EssentialReaction
truth

probabilities <- logfit.2 %>% predict(test, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
predicted.classes

mean(predicted.classes == test$EssentialReaction)

###########
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)

svm_Linear <- train(EssentialReaction ~  BetweennessCentrality+CascadeNum, data = train, method = "svmLinear",
                    trControl=trctrl,
                    preProcess = c("center", "scale"),
                    tuneLength = 10)
svm_Linear

probabilities <- svm_Linear %>% predict(test, type = "prob")
predicted.classes3 <- ifelse(probabilities > 0.5, 1, 0)
predicted.classes3

mean(predicited.class3 == test$EssentialReaction)

###########

#Now I will present an additional model that will look at Random Forest
library(randomForest)
summary(train)
rf.fit=randomForest(EssentialReaction ~  BetweennessCentrality+CascadeNum, data=train)

probabilities <- rf.fit %>% predict(test, type = "response")
predicted.classes4 <- ifelse(probabilities > 0.5, 1, 0)
predicted.classes4

mean(predicted.classes4 == test$EssentialReaction)








