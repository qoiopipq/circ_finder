import pandas as pd
import numpy as np
from pybedtools import BedTool
import sys 

#input files:
circ=sys.argv[1]#circRNAs_reads_symbols_test2.txt
junc_file=sys.argv[2]#back_spliced_junction.bed

#output files:
out_file=sys.argv[3]#junc_filtered.bed
out_filtered_df=sys.argv[4]#junc_filtered_readcounts.csv



def filter_junctions(circ, junc_file, out_file, out_filtered_df):
    df=pd.read_csv(circ, sep="\t|/", header=None, engine='python')
    df.columns=['chrom','start','end','fusion_junc','split_reads','chimeric_reads','junc_strand','gene_symbol','gene_strand','start_anno','end_anno']
    df['genomic_length']=df['end']-df['start']

    #same junction mapped to multiple genes
    fil_junc=df[(df.groupby('fusion_junc')['fusion_junc'].transform('size')==1)|
               (df.groupby('fusion_junc')['fusion_junc'].transform('size')>=2)&(df.start_anno!='Novel')]
    #exclude genes like Xist
    fil_junc=fil_junc[fil_junc.gene_symbol!='Xist']
    #only keep those genes have multiple junctions that matched with annotated splice sites
    fil_junc=fil_junc[(fil_junc.groupby('gene_symbol')['gene_symbol'].transform('size')==1)|
                (fil_junc.groupby('gene_symbol')['gene_symbol'].transform('size')>=2)&(fil_junc.start_anno!='Novel')]
    
    #Dataframe with at least 2 split reads
    fil_junc=fil_junc[fil_junc.split_reads>=2]#then used for linking fusion junction id with gene symbol
    fil_junc.to_csv(out_filtered_df)#output it for seed matching step
    filtered_junc=fil_junc['fusion_junc'].tolist()
    
    with open(out_file,'w') as junc_cand:
        with open(junc_file) as fh:
            for line in fh:
                line=line.rstrip().split('\t')
                if int(line[2])-int(line[1])<=100000: 
                    junc_id=line[3].split('/')[0]
                    if junc_id not in filtered_junc: continue
                    else: 
                        #junc_cand+=('\t'.join(line)+'\n')
                        print('\t'.join(line), file=junc_cand)
    