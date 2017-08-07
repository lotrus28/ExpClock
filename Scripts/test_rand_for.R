library(randomForest)

set.seed(1238)

sep_train_and_test <- function(df_tissue){
  ages = seq(25, 95, by=10)
  train_indices = c()
  test_indices = c()
  for (age in ages){
    # List of column indices with this age
    temp = colnames(df_tissue)[which(df_tissue['Age',] == age)]
    
    if (length(temp) == 1) {
      train_indices = c (train_indices, temp)
    } else {
      train_indices = c(train_indices, temp[1:(length(temp)*0.8)])
      test_indices = c(test_indices, temp[((length(temp)*0.8)+1):length(temp)])
    }
  }
  train_indices = train_indices[!is.na(train_indices)]
  test_indices = test_indices[!is.na(test_indices)]
  return(list(train = train_indices, test = test_indices))
}


args = commandArgs(trailingOnly=TRUE)
print(args)
p_samples = args[1]
p_metadata = args[2]
p_out = args[3]
if (substr(p_out,nchar(p_out),nchar(p_out)) != '/'){
  p_out = paste0(p_out,'/')
}

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

# Turn granular age labels into actual numbers
subj_data['Age'] = as.numeric(paste0(substr(subj_data[['Age']], 1, 1),'5'))

for (tissue in rownames(sign_feats)) {
  print(tissue)
  if (is.na(sign_feats[tissue,1]) & is.na(sign_feats[tissue,2])) next
  conf = sign_feats[tissue,1]
  tent = sign_feats[tissue,2]
  c_formula = paste0('Age ~ ', gsub(', ',' + ', conf))
  t_formula = paste0('Age ~ ', gsub(', ',' + ', conf), ' + ', gsub(', ',' + ', tent))
  
  df = read.table(paste0(p_samples, 'samples_from_',tissue,'.txt'), header = 1, row.names = 1, stringsAsFactors = F, check.names=F)
  df['Description'] = NULL
  indices = sep_train_and_test(df)
  df = df[c(strsplit(conf, ', ')[[1]], strsplit(tent, ', ')[[1]]),indices[['test']]]
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
