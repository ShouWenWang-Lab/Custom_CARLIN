input_dir=pwd;
output_dir=input_dir+"/output";
cd ../.. % Go back to the main directory
SampleList="test_Tigre_Carlin";

cfg_type="BulkDNA_Tigre";
template='Tigre';
my_CARLIN_pipeline(SampleList,cfg_type,input_dir,output_dir,template,'read_cutoff_override',3, 'read_cutoff_floor',10)
% output_all_from_summary(SampleList,output_dir,template)
%csv_reports(SampleList,output_dir,template)

SampleList="test_Tigre_Carlin,test_Tigre_Carlin";
make_allele_bank(SampleList,output_dir,template)

cd(input_dir)
