% Author: Shou-Wen Wang

This pipeline has been modified to accommodate the new sequencing design in Tigre CARLIN

To run this pipeline from Tigre carlin seqeunce, please go to the folder: @CARLIN_def
Then, manually make a copy the file CARLIN_def_Tigre.m, and rename it as CARLIN_def.m (yes, replace the current CARLIN_def.m)

In CARLIN_def_Tigre.m, I have made several changes

1, I reversed the sequence (I don't remeber why, but Li Reversed R1 and R2, and the UMI now sits at the other end of the sequence. Also, there is not all pam seqeunces are the same). I not only flipped the sequence, but also flipped the scores along the sequence, needed for alignment

2, I set obj.match_score.SecondarySequence = 20; (orignal: 30; I found that 22 is the critical number for our data)

3, I created "BulkDNA_Tigre.json" and ""BulkRNA_Tigre.json" files, which specifies relevant parameters for the sample. 
Particularly, I added a variable: QC:15.  Here I lowered the QC threshold for selecting valid UMI, as our sample may have a very low QC score. Of course, if the sample are very good, you can change the QC:15 back to QC:30.  Also, I change the pipeline such that if this QC threshold is not specified in the json file, it is assumed that it would be 30, the default setting for the cCARLIN paper. 

4, un-related to this pipeline, when I analyzed Li's data, I noticed systematic error (T->N; C->N etc) in the UMI region. I correct the error by running a simple 


To summarize, in order to run this pipeline, you can follow these steps:

Step 1: run install_CARLIN in the matlab console at the folder CARLIN_pipeline
Step 2: manully select the correct CARLIN_def.m in the folder @CARLIN_def (as discussed above)
Step 3: choose the correct sample_type file, e.g., "BulkDNA_Tigre" or "BulkRNA_Tigre", located at cfg folder. I have not yet tested other sample types. 


As a test for this pipeline, please go to the folder 'my_data/Tigre_Carlin_test', and run test_Tigre_CARLIN_pipeline.m.
You might have to change the dir_name to point to the directory where the data exists. 


