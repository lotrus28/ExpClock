library(glmnet)
set.seed(19875)

# take significant features
sign.feats = unique(strsplit(paste(na.omit(read.table('./tested_sign_features(transcript).txt', 
                        sep ='\t', header = 1, row.names = 1, stringsAsFactors = F)$Confirmed),
                        collapse = ', '),split=', ')[[1]])

# make a db [samples]x[feats/age/tissue(?)]
samp.annot = read.table('./subj_sample_annot.txt', stringsAsFactors = F, sep = '\t', header = 1, row.names = 1)
samp.annot = samp.annot[substr(rownames(samp.annot),1,1) == 'G',]
samp.annot$Age = as.numeric(paste0(substr(samp.annot[['Age']], 1, 1),'5'))
samp.annot$Tissue = gsub(' ', '', samp.annot$Tissue)

df = as.data.frame(matrix(ncol = length(sign.feats) + 3, nrow = nrow(samp.annot)))
colnames(df) = c('Age', 'Tissue', 'Sex', sign.feats)
rownames(df) = rownames(samp.annot)

df['Sex'] = samp.annot['Sex']
df['Age'] = samp.annot['Age']
df['Tissue'] = samp.annot['Tissue']

# cut data for tweaking purposes
# df = df[1:5000, 1:500]

for (t in unique(samp.annot$Tissue)){
  print(t)
  print(which(unique(samp.annot$Tissue) == t))
  tis.samp = read.table(paste0('./v7/tissue_samples/samples_from_',t,'.txt'),
                        sep = '\t', header = 1, row.names = 1, stringsAsFactors = F,
                        check.names=F)
  tis.samp$gene_id = NULL
  tis.samp = as.data.frame(t(tis.samp[rownames(tis.samp) %in% sign.feats,]))
  
  df[rownames(tis.samp), colnames(tis.samp)] <- tis.samp
  
  na_count <- apply(df[df['Tissue'] == t,], 1, function(x) sum(is.na(x)))
  df = df[!rownames(df) %in% names(na_count[na_count>0]), ]
  
}


tis.presented = as.data.frame(table(df$Tissue))
tis.presented = as.character(tis.presented[tis.presented$Freq > 200,]$Var1)
# tis.presented = c(
#   "Adipose-Subcutaneous",
#   "Adipose-Visceral(Omentum)",
#   "Artery-Aorta",
#   "Artery-Tibial",
#   "Breast-MammaryTissue",
#   "Cells-Transformedfibroblasts",
#   "Colon-Sigmoid",
#   "Colon-Transverse",
#   "Esophagus-GastroesophagealJunction",
#   "Esophagus-Mucosa",
#   "Esophagus-Muscularis",
#   "Heart-AtrialAppendage",
#   "Heart-LeftVentricle",
#   "Lung",
#   "Muscle-Skeletal",
#   "Nerve-Tibial",
#   "Pancreas",
#   "Skin-NotSunExposed(Suprapubic)",
#   "Skin-SunExposed(Lowerleg)",
#   "Stomach",
#   "Testis",
#   "Thyroid",
#   "WholeBlood"
# )
df = df[df$Tissue %in% tis.presented,]


# Even out age groups
n = 500
ages = unique(df$Age)
df_downsamp = df[0,]
for (a in ages){
  ids = rownames(df[df$Age == a,])[1:min(n,length(df[df$Age == a,2]))]
  df_downsamp = rbind(df_downsamp,df[((df$Age == a)&(rownames(df) %in% ids)),])
}


# write.table(df,'./v7/subset_v7.txt', sep = '\t', quote = F)
rm(t,i,j,na_count)



###################################################


n = nrow(df_downsamp)
m = ncol(df_downsamp)
train_rows <- sample(1:n, .66*n)
df.train <- df_downsamp[train_rows, ]
df.test <- df_downsamp[-train_rows, ]

fit.lasso <- glmnet(as.matrix(df.train[,-c(1,2,3)]), df.train$Age, family="gaussian", alpha=1)
fit.ridge <- glmnet(as.matrix(df.train[,-c(1,2,3)]), df.train$Age, family="gaussian", alpha=0)
fit.elnet <- glmnet(as.matrix(df.train[,-c(1,2,3)]), df.train$Age, family="gaussian", alpha=.5)

# fit.lasso.cv <- cv.glmnet(as.matrix(df.train[,-c(1,2,3)]), df.train$Age, type.measure="mse", family="gaussian", alpha=1)
# fit.ridge.cv <- cv.glmnet(as.matrix(df.train[,-c(1,2,3)]), df.train$Age, type.measure="mse", family="gaussian", alpha=0)
# fit.elnet.cv <- cv.glmnet(as.matrix(df.train[,-c(1,2,3)]), df.train$Age, type.measure="mse", family="gaussian", alpha=0.5)

for (i in 0:10) {
  assign(paste("fit", i, sep=""), cv.glmnet(as.matrix(df.train[,-c(1,2,3)]), df.train$Age, type.measure="mse", 
                                            alpha=i/10,family="gaussian"))
}

par(mfrow=c(3,2))

plot(fit.lasso, xvar="lambda")
plot(fit10, main="LASSO")

plot(fit.ridge, xvar="lambda")
plot(fit0, main="Ridge")

plot(fit.elnet, xvar="lambda")
plot(fit5, main="Elastic Net")

yhat0 <- predict(fit0, s=fit0$lambda.1se, newx=as.matrix(df.test[,-c(1,2,3)]))
yhat1 <- predict(fit1, s=fit1$lambda.1se, newx=as.matrix(df.test[,-c(1,2,3)]))
yhat2 <- predict(fit2, s=fit2$lambda.1se, newx=as.matrix(df.test[,-c(1,2,3)]))
yhat3 <- predict(fit3, s=fit3$lambda.1se, newx=as.matrix(df.test[,-c(1,2,3)]))
yhat4 <- predict(fit4, s=fit4$lambda.1se, newx=as.matrix(df.test[,-c(1,2,3)]))
yhat5 <- predict(fit5, s=fit5$lambda.1se, newx=as.matrix(df.test[,-c(1,2,3)]))
yhat6 <- predict(fit6, s=fit6$lambda.1se, newx=as.matrix(df.test[,-c(1,2,3)]))
yhat7 <- predict(fit7, s=fit7$lambda.1se, newx=as.matrix(df.test[,-c(1,2,3)]))
yhat8 <- predict(fit8, s=fit8$lambda.1se, newx=as.matrix(df.test[,-c(1,2,3)]))
yhat9 <- predict(fit9, s=fit9$lambda.1se, newx=as.matrix(df.test[,-c(1,2,3)]))
yhat10 <- predict(fit10, s=fit10$lambda.1se, newx=as.matrix(df.test[,-c(1,2,3)]))

mse0 <- mean((df.test$Age - yhat0)^2)
mse1 <- mean((df.test$Age - yhat1)^2)
mse2 <- mean((df.test$Age - yhat2)^2)
mse3 <- mean((df.test$Age - yhat3)^2)
mse4 <- mean((df.test$Age - yhat4)^2)
mse5 <- mean((df.test$Age - yhat5)^2)
mse6 <- mean((df.test$Age - yhat6)^2)
mse7 <- mean((df.test$Age - yhat7)^2)
mse8 <- mean((df.test$Age - yhat8)^2)
mse9 <- mean((df.test$Age - yhat9)^2)
mse10 <- mean((df.test$Age - yhat10)^2)

# nuzhen downsampling!

beta = fit.elnet$beta[, which(fit6$lambda == fit6$lambda.1se)]
beta = beta[beta != 0]
