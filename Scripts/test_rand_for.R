library(randomForest)

set.seed(1238)

args = commandArgs(trailingOnly=TRUE)
print(args)
p_samples = args[1]
p_metadata = args[2]
p_out = args[3]

# p_samples = './tissue_samples/'
# p_metadata = './subj_sample_annot.txt'
# p_out = './

sign_feats = read.table('significant_features.txt', sep = '\t', stringsAsFactors = F,
                        col.names = c('Tissue', 'Confirmed', 'Tentative'), row.names = 1)
sign_feats['RMSE Conf'] = NA
sign_feats['RSq Conf'] = NA
sign_feats['RMSE Tent'] = NA
sign_feats['RSq Tent'] = NA

subj_data = read.table(p_metadata, sep = '\t', header = 1, row.names = 1, stringsAsFactors = F)

# Turn granular age labels into actual nubers
subj_data[subj_data == '20-29'] = 25
subj_data[subj_data == '30-39'] = 35
subj_data[subj_data == '40-49'] = 45
subj_data[subj_data == '50-59'] = 55
subj_data[subj_data == '60-69'] = 65
subj_data[subj_data == '70-79'] = 75
subj_data[subj_data == '80-89'] = 85
subj_data[subj_data == '90-99'] = 95
subj_data['Age'] = as.numeric(subj_data[['Age']])

for (tissue in rownames(sign_feats)) {
  print(tissue)
  if (is.na(sign_feats[tissue,1]) & is.na(sign_feats[tissue,2])) next
  conf = sign_feats[tissue,1]
  tent = sign_feats[tissue,2]
  c_formula = paste0('Age ~ ', gsub(', ',' + ', conf))
  t_formula = paste0('Age ~ ', gsub(', ',' + ', conf), ' + ', gsub(', ',' + ', tent))
  
  df = read.table(paste0(p_samples, 'samples_from_',tissue,'.txt'), header = 1, row.names = 1, stringsAsFactors = F)
  colnames(df) = gsub('\\.','-',colnames(df))
  df = df[c(strsplit(conf, ', ')[[1]], strsplit(tent, ', ')[[1]]),-1]
  df['Age',] = subj_data[colnames(df),'Age']
  
  df = data.frame(t(df))
  # df = transform(df, class=as.numeric(as.character(df)))
  
  if (!is.na(sign_feats[tissue,1])) {
    df_to_for = df[,c(strsplit(conf, ', ')[[1]], 'Age')]
    z = randomForest(Age ~ ., data = df_to_for, ntrees = 500)
    sign_feats[tissue, 'RSq Conf'] = mean(z$rsq)
    sign_feats[tissue, 'RMSE Conf'] = mean(z$mse)**0.5
  }
  if (!is.na(sign_feats[tissue,2])) {
    z = randomForest(Age ~ ., data = df, ntrees = 500)
    sign_feats[tissue, 'RSq Tent'] = mean(z$rsq)
    sign_feats[tissue, 'RMSE Tent'] = mean(z$mse)**0.5
  }
}

write.table(sign_feats, file = paste0(p_out,'tested_sign_features.txt'), quote = F, sep = '\t')
