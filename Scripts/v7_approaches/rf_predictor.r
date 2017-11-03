# compare speed of aging between diff-t tissues and sexes

# Total: 18 tissues

# Sexes:
# F
# M

# Create df:
# tis | sex | age | pred

# Read:
# Subject-Sample annotation -- fill (Subject | Age)

# Read:
# Confident tissue predictors

# Read:
# Tissue1_data -- use predictors to create RF, predict age, fill (| Age_Tissue1|) column
library(randomForest)
library(ggplot2)
library(gridExtra)
library(lme4)
library(arm)
require(reshape)

# tis = c('Kidney-Cortex')
tis = c(
    "Adipose-Subcutaneous",
    "Adipose-Visceral(Omentum)",
    "Artery-Aorta",
    "Artery-Tibial",
    "Breast-MammaryTissue",
    "Cells-Transformedfibroblasts",
    "Colon-Sigmoid",
    "Colon-Transverse",
    "Esophagus-GastroesophagealJunction",
    "Esophagus-Mucosa",
    "Esophagus-Muscularis",
    "Heart-AtrialAppendage",
    "Heart-LeftVentricle",
    "Lung",
    "Muscle-Skeletal",
    "Nerve-Tibial",
    "Pancreas",
    "Skin-NotSunExposed(Suprapubic)",
    "Skin-SunExposed(Lowerleg)",
    "Stomach",
    "Testis",
    "Thyroid",
    "WholeBlood"
)

sex = c('F', 'M')
conf = read.table('./v7/tested_sign_features(transcript).txt', stringsAsFactors = F, sep ='\t', header = 1, row.names = 1)

subj_samp = read.table('./v7/subj_sample_annot.txt', stringsAsFactors = F, sep = '\t', header = 1, row.names = 1)
subj_samp = subj_samp[substr(rownames(subj_samp),1,1) == 'G',]
subj_samp$Subject.ID = regmatches(rownames(subj_samp),regexpr("^GTEX-\\w+",rownames(subj_samp)))
subj_samp$Prediction = NA
subj_samp$Age = as.numeric(paste0(substr(subj_samp[['Age']], 1, 1),'5'))

subjects = as.data.frame(matrix(ncol = 2, nrow = length(unique(subj_samp$Subject.ID))))
colnames(subjects) = c('Age','ID')
rownames(subjects) = unique(subj_samp$Subject.ID)
for (i in rownames(subjects)){
  subjects[i,1] = subj_samp[subj_samp$Subject.ID == i,'Age'][1]
  subjects[i,2] = i
}
  
# Downsample
# 53/52/115/235/226/21
n = 40
ages = unique(subj_samp$Age)
downsamp = subj_samp[0,]
for (a in ages){
  ids = subjects[subjects$Age == a,2][1:min(n,length(subjects[subjects$Age == a,2]))]
  downsamp = rbind(downsamp,subj_samp[((subj_samp$Age == a)&(subj_samp$Subject.ID %in% ids)),])
}

subj_samp = downsamp

for (t in tis) {
  set.seed(1238)
  
  tis_samp = read.table(paste0('./v7/tissue_samples/samples_from_', t,'.txt'), stringsAsFactors = F, 
                        check.names=F, header = 1, row.names = 1)
  
  tis_samp$gene_id <- NULL
  tis_samp = as.data.frame(t(tis_samp))
  tis_samp$Age = subj_samp[rownames(tis_samp),'Age']
  conf_tis = strsplit(conf[t,1], split = ', ')[[1]]
  if (is.na(conf_tis[1])){
    next
  }
  
  tis_samp = tis_samp[,c(conf_tis,'Age')]
  tis_samp = na.omit(tis_samp)
  
  rf = randomForest(Age ~ ., data = tis_samp, ntrees = 128)
  rmse = (sum(as.vector(rf$predicted - tis_samp$Age)**2)/length(tis_samp$Age))**0.5
  print(paste0('RMSE for ', t, ' is ', round(rmse,2), 'yrs'))
  
  subj_samp[names(rf$predicted),'Prediction'] = rf$predicted
}

subj_samp = subj_samp[(!(is.na(subj_samp$Prediction))),]

subj_preds = as.data.frame(matrix(ncol = length(c(unique(subj_samp$Tissue))) + 2, nrow = length(unique(subj_samp$Subject.ID))))
colnames(subj_preds) = c(unique(subj_samp$Tissue), 'Age', 'Prediction')
rownames(subj_preds) = unique(subj_samp$Subject.ID)

for (i in rownames(subj_preds)){
  temp = subj_samp[subj_samp$Subject.ID == i,]
  for (j in rownames(temp)){
    subj_preds[i,temp[j,'Tissue']] = temp[j,'Prediction']
    subj_preds[i,'Sex'] = temp[1,'Sex']
    subj_preds[i,'Age'] = temp[1,'Age']
  }
}
colnames(subj_preds) = make.names(colnames(subj_preds), unique = T)

x = rfImpute(Age ~ ., data = subj_preds[,!(names(subj_preds) %in% c('Rand_age', 'Error','Subject','Sex','Prediction'))], ntree = 128)
y = randomForest(Age ~ ., data = x, ntrees = 128)

subj_preds$Prediction = y$predicted
subj_preds$Rand_age = subj_preds$Age + runif(nrow(subj_preds), -5, 5)

MLexamp <- lmer(Age ~ Prediction + Tissue + (1 + Prediction | Sex), data = subj_samp)
display(MLexamp)
coef(MLexamp)
############
############
############
############
############
############
df = subj_preds
sc_plot <- ggplot(df, aes(df$Prediction, df$Rand_age, color = df$Sex)) +
  geom_point(alpha=.3, size = 3) +
  scale_color_manual(values = c('#109910', '#991010')) +
  theme_classic() +
  theme(legend.position='right') +
  labs(x = 'Prediction', y = 'Age')
sc_plot

# subj_samp$Error = subj_samp$Age - subj_samp$Prediction
# boxplot(subj_samp[subj_samp['Sex'] == 'F','Error'],subj_samp[subj_samp['Sex'] == 'M','Error'])
# subj_samp$Rand_age = subj_samp$Age + runif(nrow(subj_samp), -5, 5)

############
# + Subject.ID

# MLexamp <- glm(Age ~  Prediction + Tissue*Sex, data = subj_samp)
# x = (as.data.frame(MLexamp$coefficients))

############

temp = subj_preds[,!(names(subj_preds) %in% c('Subject','Sex','Rand_age','Error'))]
temp = temp - subj_preds$Age
boxplot(temp)


tis_cor = list()
for (i in colnames(temp)){
  tis_cor[i] = cor(temp[!is.na(temp[i]),c(i,'Prediction')])[1,2]
}

for (i in rownames(x)){
  x[i,] = rank(x[i,])
}



p1 <- ggplot(data = melt(x[x$Age == 25,]), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable)) +
  theme(legend.position="none") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylim(20,65)
p2 <- ggplot(data = melt(x[x$Age == 35,]), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable)) + 
  theme(legend.position="none") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylim(20,65)
p3 <- ggplot(data = melt(x[x$Age == 45,]), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable)) +
  theme(legend.position="none") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylim(30,70)
p4 <- ggplot(data = melt(x[x$Age == 55,]), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable)) +
  theme(legend.position="none") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylim(30,65)
p5 <- ggplot(data = melt(x[x$Age == 65,]), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable)) +
  theme(legend.position="bottom") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p6 <- ggplot(data = melt(x[x$Age == 75,]), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable)) +
  theme(legend.position="none") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2, widths = c(6,10,6), heights = c(4, 4))
