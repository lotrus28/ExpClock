# install.packages('Boruta')
# install.packages('randomForest')
library(Boruta)
library(randomForest)

args = commandArgs(trailingOnly=TRUE)
print(args)
p_samples = args[1]
p_subjects = args[2]
out_folder = args[3]

# p_samples = 'samples_from_Brain-Hippocampus.txt'
# p_subjects = 'subj_sample_annot.txt'

# Read RPKM data for a tissue
# Read subject metadata: age and sex
df = read.table(p_samples, sep = '\t', header = 1, row.names = 1, stringsAsFactors = F)
subj_data = read.table(p_subjects, sep = '\t', header = 1, row.names = 1, stringsAsFactors = F)

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

# Replace dots with dashes in column names
colnames(df) = gsub('\\.','-',colnames(df))

# Add sex and age columns to sample table
df['Sex',] = 0
df['Age',] = 0
for (sample in colnames(df)[2:ncol(df)]) {
  df['Sex',sample] = subj_data[[sample,'Sex']]
  df['Age',sample] = subj_data[[sample,'Age']]
}

# Separate training set
train_df = df[,1:(ncol(df)*0.8)]
test_df = df[,((ncol(df)*0.8)+1):ncol(df)]

# Remove gene description field from Boruta input
temp = train_df[, 2:ncol(train_df)]
# Erase sex variable
# temp = temp[!(rownames(temp) == 'Sex'),]
# Take only 1'000 features
temp = temp[c(rownames(temp)[rownames(temp) != 'Age'][1:1000],'Age'),]

# Remove all zero rows
temp = temp[!(apply(temp, 1, function(y) all(y == 0))),]

# Transpose so that features are in columns
# and samples are in rows
df_to_bor = t(temp)
class(df_to_bor) = 'numeric'
df_to_bor = as.data.frame(df_to_bor)

# Perform Boruta random forest feature selection
set.seed(124)
# ntree is passed to randomForest later
bor_brain = Boruta(Age ~ ., data = df_to_bor, doTrace = 2, ntree = 300, maxRuns = 800)

# View(as.data.frame(bor_brain$finalDecision))
# getNonRejectedFormula(bor_brain)
# getConfirmedFormula(bor_brain)
# getSelectedAttributes(bor_brain)
# View(attStats(bor_brain))
# x = randomForest(df_to_bor[,tent],df_to_bor$Age)

tent = names(bor_brain$finalDecision[bor_brain$finalDecision == 'Tentative'])
conf = names(bor_brain$finalDecision[bor_brain$finalDecision == 'Confirmed'])
rej  = names(bor_brain$finalDecision[bor_brain$finalDecision == 'Rejected'])

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