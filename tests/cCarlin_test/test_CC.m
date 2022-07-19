% input_dir=pwd;
% output_dir=input_dir+"/output";
% cd ../.. % Go back to the main directory
% SampleList="CC-DNA-test";
% 
% cfg_type="BulkRNA_12UMI";
% template='cCARLIN';
% my_CARLIN_pipeline(SampleList,cfg_type,input_dir,output_dir,template,'read_cutoff_override',1, 'read_cutoff_floor',1)
% %output_selected_from_summary(SampleList,output_dir,template)
% %csv_reports(SampleList,output_dir,template)
% 
% % SampleList="CC-DNA-test";
% %merge_samples(SampleList,output_dir,template)
% %make_allele_bank(SampleList,output_dir,template)
% 
% cd(input_dir)



% input_dir=pwd;
% output_dir=input_dir+"/output";
% SampleList="CC_RNA_poor";
% cd ../.. % Go back to the main directory
% cfg_type="BulkRNA_12UMI";
% template='cCARLIN';
% my_CARLIN_pipeline(SampleList,cfg_type,input_dir,output_dir,template,'read_cutoff_override',1, 'read_cutoff_floor',1)
% cd(input_dir)


% input_dir=pwd;
% output_dir=input_dir+"/output";
% SampleList="CC_DNA_poor";
% cd ../.. % Go back to the main directory
% cfg_type="BulkDNA_12UMI";
% template='cCARLIN';
% my_CARLIN_pipeline(SampleList,cfg_type,input_dir,output_dir,template,'read_cutoff_override',1, 'read_cutoff_floor',1)
% cd(input_dir)


%% single-cell test
input_dir=pwd;
output_dir=input_dir+"/output";
SampleList="scLimeCat_test";
cd ../.. % Go back to the main directory
cfg_type="scLimeCat";
template='cCARLIN';
my_CARLIN_pipeline(SampleList,cfg_type,input_dir,output_dir,template,'read_cutoff_override',1, 'read_cutoff_floor',1)
cd(input_dir)