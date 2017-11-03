import pandas as pd

def main():
    subj_attr = '../v7/GTEx_v7_Annotations_SubjectPhenotypesDS.txt'
    samp_attr = '../v7/GTEx_v7_Annotations_SampleAttributesDS.txt'

    samples = pd.read_table(samp_attr, sep='\t', header=0, index_col=None)
    subj = pd.read_table(subj_attr, sep='\t', header=0, index_col=None)

    # end_table - sample ID, tissue, subject, age, sex
    columns = ['Sample ID', 'Tissue', 'Age', 'Sex']
    end_table = pd.DataFrame(columns=columns)

    for i in samples.index:
        try:
            new_line = pd.DataFrame(data=[samples.ix[i, 'SAMPID'], samples.ix[i, 'SMTSD'], 0, 0], index=columns).T
            temp = subj[subj['SUBJID'] == '-'.join(new_line.ix[new_line.index[0], 'Sample ID'].split('-')[:2])]
            new_line['Age'] = temp.ix[temp.index[0], 'AGE']
            try:
                # v6p
                new_line['Sex'] = str(temp.ix[temp.index[0], 'GENDER']).replace('1', 'M').replace('2', 'F')
            except KeyError:
                # v7
                new_line['Sex'] = str(temp.ix[temp.index[0], 'SEX']).replace('1', 'M').replace('2', 'F')
            end_table = end_table.append(new_line, ignore_index=True)
        except:
            print('Over')
            continue

    end_table.to_csv('subj_sample_annot.txt', sep='\t', index=None)

if __name__ == '__main__':
    main()

