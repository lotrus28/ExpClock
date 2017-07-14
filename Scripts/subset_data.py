import pandas as pd
import sys

# 'GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct'
global rpkm_path
rpkm_path = sys.argv[1]
global out_path
out_path = sys.argv[2]

if not(out_path[-1] == '/'):
    out_path += '/'

def truncate_rpkm_table(row_start, row_end, col_start, col_end, rpkm_path):

    trunc_data = 'truncated_rpkm_{}-{}x{}-{}.txt'.format(str(row_start), str(row_end),str(col_start), str(col_end))
    with open(trunc_data, 'a') as wr:
        with open(rpkm_path) as f:
            i = 0
            for line in f:

                if i < (row_start-1):
                    i += 1
                    continue
                if i == row_end:
                    break
                new_row = '\t'.join(line.split('\t')[col_start:(col_end+1)]) + '\n'
                wr.write(new_row)
                i += 1

def extract_tissue_samples(tissue):
    annot = pd.read_table('subj_sample_annot.txt', sep='\t', header=0, index_col=None)
    annot_trunc = annot[annot['Tissue'] == tissue]
    trunc_data = '{}samples_from_{}.txt'.format(out_path, tissue.replace(' ',''))

    with open(trunc_data, 'a') as wr:

        with open(rpkm_path) as f:

            i = 0
            for line in f:
                if i == 2:
                    header = line.split('\t')
                    start = header[:2]
                    end = [x for x in header if x in list(annot_trunc['Sample ID'])]
                    columns = [header.index(x) for x in header if x in list(annot_trunc['Sample ID'])]
                    header = '\t'.join(start + end) + '\n'
                    wr.write(header)
                if i > 2:
                    new_row = line.split('\t')
                    row_to_write = [new_row[0],new_row[1]]
                    for k in columns:
                        row_to_write.append(new_row[k])
                    row_to_write = '\t'.join(row_to_write) + '\n'
                    wr.write(row_to_write)
                i += 1

def main():
    # rpkm = pd.read_table(rpkm_path, engine='python',  nrows=20)
    # truncate_rpkm_table(3, 1000, 0, 1000, rpkm_path)

    #tissue = 'Brain - Hippocampus'
    df = pd.read_table('subj_sample_annot.txt', header = 0)
    tissues = list(set(df['Tissue']))
    for tissue in tissues:
        extract_tissue_samples(tissue)
        print(tissue + ' done')

if __name__ == '__main__':
    main()
