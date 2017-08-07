library(ggplot2)
library(gplots)
library(d3heatmap)
library(boot)

# Draw an html heatmap of sig.corrs
# in: paths to the tables of pvals and corrs
# out: d3heatmap (default)
# Might produce a heatmap.2
draw_matrix_plot <- function(p_pval = NA ,p_corr =NA, matr = NA, h2 = NULL){
  
  if (!(is.na(p_pval))) {
    sign.cor <- get_sig_cor(p_pval, p_corr)
  } else {
    sign.cor <- matr
  }
  
  my_palette <- colorRampPalette(c("red", "white", "green"))(n = 299)
  
  if(!is.null(h2)){
    sign.cor [] <- round(sign.cor[]*100,0)
    heatmap.2(sign.cor,
              cellnote = sign.cor,  # same data set for cell labels
              main = "Correlation", # heat map title
              notecol="black",      # change font color of cell labels to black
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(12,9),     # widens margins around plot
              col=my_palette,       # use on color palette defined earlier
              dendrogram="none",     # only draw a row dendrogram
              Rowv=TRUE, 
              Colv=TRUE,
              key=FALSE
    )
  } else {
    a <- d3heatmap(cor.table,
                   col=my_palette,
                   dendrogram="both",
                   Rowv=TRUE, 
                   Colv=TRUE,
                   width = 7000,
                   height = 3000,
                   digits = 1,
                   show_grid = F
    )
    b <- d3heatmap(cor.table,
                   col=my_palette,
                   dendrogram="none",
                   width = 7000,
                   height = 3000,
                   digits = 1,
                   show_grid = F
    )
    b$x$matrix <- a$x$matrix
    b$x$image  <- a$x$image
    b   
  }
}

pval_boot <- function(table, num_boot){
  true_cor = cor(table)[1,2]
  
  pear_cor <- function(data,ind){
    data = data[ind,]
    return(cor(data, method = 'pearson')[1,2])
  }
  
  boots <- boot(data = table, R=num_boot, statistic=pear_cor)$t[,1]
  
  pV = length(boots[abs(boots) > abs(true_cor)]) / length(boots)
  
  return(pV)
}

############################################

args = c('./tissue_samples/', './tested_sign_features.txt')

tissue_folder = args[1]
sign.genes_path = args[2]

tissue_samples = list.files(path = tissue_folder)
sign.genes = read.csv(sign.genes_path, sep = '\t', row.names = 1, stringsAsFactors = F)

conf = sign.genes[['Confirmed']]
conf = conf[!is.na(conf)]
all_genes = unique(strsplit(toString(conf), split = ', ')[[1]])

cor.table = data.frame(matrix(NA, ncol = length(tissue_samples), nrow = length(all_genes)))
colnames(cor.table) = rownames(sign.genes)
rownames(cor.table) = all_genes

subj_samp = read.table('subj_sample_annot.txt', sep = '\t', header = 1, row.names = 1)
subj_samp['Age'] = as.numeric(paste0(substr(subj_samp[['Age']], 1, 1),'5'))

if(!file.exists('cor_table.txt')){
  
  for (tis in colnames(cor.table)){
    sample = read.csv(file = paste0(tissue_folder, 'samples_from_', tis, '.txt'), sep = '\t', row.names = 1, header = 1 )
    sample = sample[c(all_genes),]
    colnames(sample) = gsub('\\.', '-', colnames(sample))
    for (g in rownames(cor.table)){
      temp = sample[g,]
      desc = temp[g,'Description']
      temp$Description <- NULL
      temp['Age',] = subj_samp[colnames(temp), 'Age']
      temp = t(temp)
      cor.table[g,tis] = cor(temp, method = 'pearson')[1,2]
    }
  }
  
  cor.table = cor.table[,!apply(is.na(cor.table), 2, all)]
  
  write.table(cor.table, file = 'cor_table.txt', quote = F, sep = '\t', col.names=NA)
} else {
  cor.table = read.table('cor_table.txt', sep = '\t', stringsAsFactors = F, header = 1, row.names = 1, check.names=F)
}

cor.table = cor.table[colSums(!is.na(cor.table)) > 0]
#############

pv_table = cor.table
pv_table = pv_table[colSums(!is.na(pv_table)) > 0]
pv_table[] = NA


for (tis in colnames(cor.table)) {
  
  sample = read.csv(file = paste0(tissue_folder, 'samples_from_', tis, '.txt'), sep = '\t', row.names = 1, header = 1 )
  sample = sample[c(all_genes),]
  colnames(sample) = gsub('\\.', '-', colnames(sample))
  sample$Description <- NULL
  sample['Age',] = subj_samp[colnames(sample), 'Age']
  sample = as.data.frame(t(sample))
  
  for (gen in colnames(sample)) {
    # print(gen)
    if (gen == 'Age') {
      pv_table[gen,tis] = NA
    } else {
      temp = sample[,c(gen,'Age')]
      pv_table[gen,tis] = pval_boot(temp,100)
      
    }
  }
  print(paste0('----- Completed tissue ',tis))
}

write.table(pv_table, file = 'cor_pval_table.txt', quote = F, sep = '\t', col.names=NA)

