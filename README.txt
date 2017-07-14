╔═══════════════════════════════════════════════════════════════════════════════════════════════╗
║ 																								                                              ║
║ Minor errors might be present. Pay extra attention to passed arguments and hardcoded paths!	  ║
║																								                                                ║
║ To work with scripts you need to download these files from GTEx site:							            ║
║ - GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.7z										                    ║
║ - GTEx_Data_V6_Annotations_SampleAttributesDS.txt												                      ║
║ - GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt											                      ║
║																								                                                ║
║ Data division into testing and training sets is not supported yet.							              ║
╚═══════════════════════════════════════════════════════════════════════════════════════════════╝

= get_sign_features.sh
Maaster bash script for separating data in tissue tables, finding significant features and evaluating them.
The script is not optimized for parallel computation.
You need to specify script-containing folder inside of this file

= Scripts/
Contains scripts required by `get_sign_features.sh`

= subj_sample_annot.txt
Contains  metadata available on GTEx samples: tissue, granular age and sex.
The file is produced by './Scripts/join_sample_and_subj_annot.py from sample and subject annotations.

= significant_features.txt
1st column: tissue
2nd column: Boruta confirmed features
3d column: Boruta tentative features
The file is produced by './Scripts/boruta_clock.R'

= tested_sign_features.txt
Same as significant_features.txt but has 4 columns describing prediction quality: RMSE is root mean squere error and RSq is the fraction of variance explained. Some trees include only confirmed features and others — both confirmed and tentative.
