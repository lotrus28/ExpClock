# install.packages('Boruta')
# install.packages('randomForest')
library(Boruta)

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

set.seed(124)

args = commandArgs(trailingOnly=TRUE)
print(args)
p_samples = args[1]
p_subjects = args[2]
out_folder = args[3]
if (substr(out_folder,nchar(out_folder),nchar(out_folder)) != '/'){
  out_folder = paste0(out_folder,'/')
}

# p_samples = './tissue_samples/samples_from_Bladder.txt'
# p_subjects = 'subj_sample_annot.txt'

# Read RPKM data for a tissue
# Read subject metadata: age and sex
# Using `row.names = NULL` due to some strange bug in R CMD BATCH
df = read.table(p_samples, sep = '\t', header = 1, row.names = NULL, stringsAsFactors = F, check.names=F)
rownames(df) = df[['Name']]
df['Name'] = NULL
df$Description <- NULL
subj_data = read.table(p_subjects, sep = '\t', header = 1, row.names = 1, stringsAsFactors = F)

# Turn granular age labels into actual numbers
subj_data['Age'] = as.numeric(paste0(substr(subj_data[['Age']], 1, 1),'5'))

# Add sex and age columns to sample table
# df['Sex',] = 0
df['Age',] = 0
for (sample in colnames(df)) {
  # df['Sex',sample] = subj_data[[sample,'Sex']]
  df['Age',sample] = subj_data[[sample,'Age']]
}

# Separate training set
indices = sep_train_and_test(df)
train_df = df[,indices[['train']]]
# test_df = df[,indices[['test']]]

# Erase sex variable
# temp = temp[!(rownames(temp) == 'Sex'),]

# Remove all zero rows
train_df = train_df[!(apply(train_df, 1, function(y) all(y == 0))),]

# Transpose so that features are in columns
# and samples are in rows
df_to_bor = t(train_df)
class(df_to_bor) = 'numeric'
df_to_bor = as.data.frame(df_to_bor)
train_ages = df_to_bor['Age']
df_to_bor['Age'] = NULL
df_to_bor = as.matrix(df_to_bor)

# df_to_bor = df_to_bor[,1:4000]

# Perform Boruta random forest feature selection
total = 1

bor_brain = Boruta(x = df_to_bor, y = train_ages[[1]], doTrace = 2, maxRuns = 3000, holdHistory = F)
tent = names(bor_brain$finalDecision[bor_brain$finalDecision == 'Tentative'])
conf = names(bor_brain$finalDecision[bor_brain$finalDecision == 'Confirmed'])
rej  = names(bor_brain$finalDecision[bor_brain$finalDecision == 'Rejected'])

df_to_bor = df_to_bor[,c(tent, conf)]

# View(as.data.frame(bor_brain$finalDecision))
# getNonRejectedFormula(bor_brain)
# getConfirmedFormula(bor_brain)
# getSelectedAttributes(bor_brain)
# View(attStats(bor_brain))
# x = randomForest(df_to_bor[,tent],df_to_bor$Age)

filename = strsplit(p_samples, split='/|_|\\.')[[1]]
tissue = filename[length(filename)-1]
if (length(rej) == (ncol(df_to_bor) - 1) ){
  print('No features significant')
  write(c(tissue, NA, NA), file = 'significant_features.txt', append = T, sep = '\t', ncolumns = 3)
} else{
  print(paste0(length(tent), ' tentative features found'))
  print(paste0(length(conf), ' confirmed features found'))
  if (length(conf) == 0) {
    conf = NA
  } else {
    conf = paste(conf, collapse = ", ")
  }
  if (length(tent) == 0) {
    tent = NA
  } else {
    tent = paste(tent, collapse = ", ")
  }
  write(c(tissue, conf, tent), file = paste0(out_folder, 'significant_features.txt'), append = T, sep = '\t', ncolumns = 3)
}