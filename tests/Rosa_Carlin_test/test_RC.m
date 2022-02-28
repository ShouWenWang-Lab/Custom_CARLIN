input_dir=pwd;
output_dir=input_dir+"/output";
cd ../.. % Go back to the main directory
SampleList="RC-DNA-test";

cfg_type="BulkDNA_Rosa_14UMI";
template='Rosa';
my_CARLIN_pipeline(SampleList,cfg_type,input_dir,output_dir,template,'read_cutoff_override',1, 'read_cutoff_floor',1)
% output_all_from_summary(SampleList,output_dir,template)
%csv_reports(SampleList,output_dir,template)

SampleList="RC-DNA-test";
%merge_samples(SampleList,output_dir,template)
make_allele_bank(SampleList,output_dir,template)

cd(input_dir)