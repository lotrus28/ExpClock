╔═══════════════════════════════════════════════════════════════════════════════════════════════╗
║
║ Minor errors might be present. Pay extra attention to passed arguments and hardcoded paths!
║
║ To work with scripts you need to download these files from GTEx site:
║ - GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.7z
║ - GTEx_Data_V6_Annotations_SampleAttributesDS.txt
║ - GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt
║
║ Data division into testing and training sets is not supported yet.
╚═══════════════════════════════════════════════════════════════════════════════════════════════╝

= subset_data.py
Separates the giant RPKM table for all samples into numerous tables containing samples form one organ.
in: RPKM table path, output folder path
out: RPRM tables for different organs

= boruta_clock.R
Using Boruta to identify significant features (genes) in RPKM tables. The script produces 'significant_features.txt'
in: path to RPKM table, path to subject-sample combined annotation, output folder
out: table of features of confirmed and tentative significance as determined by Boruta

= test_rand_for.R
Building random forest with Boruta-confirmed and -tentative features.
Output — 'tested_sign_features.txt'. It is the same as 'significant_features.txt' but also contains 4 columns describing quality of predictions: RMSE is root mean squere error and RSq is the fraction of variance explained. Some trees include only confirmed features and others — both confirmed and tentative.
in: folder with RPKM tables, subject-sample metadata table, output folder
out: table of confirmed-tentative features per tissue added with prediction quality
