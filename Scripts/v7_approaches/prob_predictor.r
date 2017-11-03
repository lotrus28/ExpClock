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

conf.all = na.omit(unique(unlist(strsplit(conf$Confirmed, split = ', '))))



tis.probs = as.data.frame(matrix(data = NA, ncol = 8, nrow = 0))
colnames(tis.probs) = c('Subject', 'Tissue', 'a25', 'a35', 'a45', 'a55', 'a65', 'a75')
for (t in tis[6:length(tis)]) {
  start = Sys.time()
  set.seed(1238)
  
  tis_samp = read.table(paste0('./v7/tissue_samples/samples_from_', t,'.txt'), stringsAsFactors = F, 
                        check.names=F, header = 1, row.names = 1)
  
  tis_samp$gene_id <- NULL
  tis_samp = as.data.frame(t(tis_samp))
  tis_samp$Age = subj_samp[rownames(tis_samp),'Age']
  
  # conf_tis = strsplit(conf[t,1], split = ', ')[[1]]
  conf_tis = conf.all
  
  if (is.na(conf_tis[1])){next}
  
  insufficient_data = F
  for (i in table(tis_samp$Age)){if (i < 5) {insufficient_data = T}}
  if (insufficient_data){next}
  
  tis_samp = tis_samp[,c(conf_tis,'Age')]
  tis_samp = na.omit(tis_samp)
  
  for (s in rownames(tis_samp)){
    ID = regmatches(s,regexpr("^GTEX-\\w+",s))
    tis.probs[nrow(tis.probs) + 1, 'Subject'] = ID
    tis.probs[nrow(tis.probs),'Tissue'] = t
    
    feats = colnames(tis_samp)[1:(length(tis_samp)-1)]
    ages = c(25,35,45,55,65,75)
    
    for (a in ages){
      age_probs = c()
      for (tra in feats){
        tra_dist = tis_samp[tis_samp$Age == a, tra]
        med_tra = median(tra_dist)
        if (tis_samp[s,tra] < med_tra){
          prob = (sum(tis_samp[s,tra] >= tra_dist) + 1) / (length(tra_dist) + 1)
        } 
        if (tis_samp[s,tra] > med_tra){
          prob = (sum(tis_samp[s,tra] <= tra_dist) + 1) / (length(tra_dist) + 1)
        }
        if (tis_samp[s,tra] == med_tra){
          prob = max(0.5, (sum(tra_dist == med_tra) / length(tra_dist)))
        }
        age_probs = c(age_probs, prob)
        tis.probs[nrow(tis.probs), paste0('a',a)] = sum(log(age_probs))
      }
    }
  }
  end = Sys.time()
}
print(end - start)

for (i in rownames(tis.probs)){
  tis.probs[i,'Age'] = subjects[tis.probs[i,'Subject'], 'Age']
  tis.probs[i,'Predicted'] = substr(
      colnames(tis.probs)[which(tis.probs[i,] == max(tis.probs[i,c(3,4,5,6,7,8)]))],
      2,3)
}
tis.probs$Predicted = as.numeric(tis.probs$Predicted)
rmse = (sum((tis.probs$Predicted - tis.probs$Age)**2)/length(tis.probs$Age))**0.5

for (i in rownames(subjects)){
  subjects[i,'Predicted'] = mean(tis.probs[tis.probs$Subject == i,'Predicted'])
}

rmse = (sum(na.omit(subjects$Predicted - subjects$Age)**2)/length(subjects$Age))**0.5

write.table(tis.probs, 'probabilities.txt', sep ='\t', quote = F, row.names = F)
