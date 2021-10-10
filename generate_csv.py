import pandas as pd

import argparse
# parse cmd line arguments
parser = argparse.ArgumentParser(description="Generate csv")
parser.add_argument(
    "--path",
    type=str,
    default=".",
    help="raw data path",
)

data_dir = parser.parse_args().path

def read_from_matlab_output(source_path):
    sample_list=open(f'{source_path}/sample_names.txt').read().split('\n')[:-1]
    variable_list=open(f'{source_path}/variable_names.txt').read().split('\n')[:-1]
    df=pd.read_csv(f'{source_path}/result_2.csv',names=variable_list)
    df['sample']=sample_list
    cols = df.columns.tolist()
    df=df[cols[-1:]+cols[:-1]]
    source=[xx[:2] for xx in sample_list]
    df['source']=source
    df.to_csv(f'{source_path}/refined_results.csv',float_format='%.2f')
    return df




df_out=read_from_matlab_output(data_dir)
