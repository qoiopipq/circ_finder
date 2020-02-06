import pandas as pd
import numpy as np
from pybedtools import BedTool
import regex as re
import sys 
import io
import os
import timeit
from Bio import SeqIO
from collections import defaultdict

#input files:
junc_cand=sys.argv[1]#filtered_junc.bed from filter_junc step
ref_exon=sys.argv[2]#mouse_exons.bed
genome=sys.argv[3]#genome.fa
fil_junc=sys.argv[4]#filtered_readcounts.csv from filter_junc step
seeds=sys.argv[5]#mmu_miRNAs_seeds.fa

#output files:
out_file=sys.argv[6]#circ_seed_matches.csv
   
def extract_exon_seq(junc_cand, ref_exon,genome, fil_junc, seeds, out_file):   
    junc_cand=BedTool(junc_cand)#turn this into a BedTool object
    exons_ref=BedTool(ref_exon)#ref exons in bed format from UCSC browser
    circ_exons=junc_cand.intersect(exons_ref)
    item_list=[]
    circ_exon_uniq=''
    for item in circ_exons:
            item_s=str(item)
            if item_s not in item_list:
                item_list.append(item_s)
                circ_exon_uniq+=(str(item)+'\n')
    circ_exon_uniq=BedTool(circ_exon_uniq,from_string=True)#turn this into a BedTool object
    with open(genome) as fasta:
        #save the intermediate file in the current dir
        circ_exons_seqs=circ_exon_uniq.sequence(fi=fasta,name=True).save_seqs('circ_cand_seq.fa')
        
        
    #Concatenate sequences for same junction
    sequence_map = defaultdict(str)
    

    for cur_record in SeqIO.parse(os.getcwd()+'/circ_cand_seq.fa', "fasta"):
        sequence_map[cur_record.name] += str(cur_record.seq)

    
    conc_circ_cand_seq=''
    for k,v in sequence_map.items():
        conc_circ_cand_seq+=('>'+k+'\n'+v+'\n')#conc_circ_cand_seq.fa
            #print('>'+k+'\n'+v', file=conc_circ_cand_seq)
    
    return (seed_match(seeds, io.StringIO(conc_circ_cand_seq), fil_junc, out_file))

def seed_match(seeds, conc_circ_cand_seq, fil_junc, out_file):
    #Link fusion junction id with gene symbol
    fil_junc=pd.read_csv(fil_junc,header=0, index_col=0)
    fil_junc=fil_junc.sort_values(by=['split_reads','chimeric_reads'],ascending=False)
    #new data frame with 7 columns:
    junc_symbol=fil_junc[['fusion_junc','gene_symbol','split_reads','chimeric_reads','chrom','start_anno','end_anno']]
    #convert columns:'fusion_junc' and 'gene_symbol' into a dictionary
    junc_symbol_convert=junc_symbol.set_index('fusion_junc').to_dict()
    junc_symbol_convert=junc_symbol_convert['gene_symbol']  
    
    df_meta=pd.DataFrame({'fusion_junc':[],'length':[],'gene_symbol':[],'num_miRNA_hits':[],
                'miRNA_list':[],'max_seed_matches':[],'max_miRNA_list':[]})
    
    for index, seq_record in enumerate(SeqIO.parse(conc_circ_cand_seq,'fasta')):
        miR_list=[]#store # of matching miRNAs
        max_match=0#for selecting local max matches
        max_miR=dict()#for selecting global max match
        for i, record in enumerate(SeqIO.parse(seeds,'fasta')):
            h=record.name
            s=record.seq
            rc=s.reverse_complement()
            #check for seed matches in circ
            seed_matches=re.findall(str(rc), str(seq_record.seq),re.I,overlapped=True)
            #show miRNAs with more than 1 seed match 
            num_matches=len(seed_matches)
            if seed_matches and num_matches<=1:
                continue
                #print(str(seq_record.name).split('/')[0])
                #print(h+', Seed: '+str(s)+', Matches =',num_matches)
            elif seed_matches and num_matches>1:
                miR_list.append(h)
                #record the max seed matches for each circRNA candidate
                if max_match<=num_matches:
                    max_match=num_matches
                    max_miR[h]=max_match
                circ_id=str(seq_record.name).split('/')[0]      
                #add in one more column, lengths of circRNA exonic sequences🌟
                circ_exon_len=len(seq_record.seq)
                gene_symbol=''.join([v for k,v in junc_symbol_convert.items() if k==circ_id])#!!!!!!!!
                max_miR_IDs=', '.join([k for k,v in max_miR.items() if v==max(max_miR.values())])
        #add this if statement to exclude incorrect seed matching 🌟
        if max_match==0: continue
        #print(circ_id,gene_symbol)
        #print('Max seed matches:', max_match,str(max_miR_IDs))
        #print('Matched miRNAs: ', len(miR_list),', '.join(miR_list))
        df_meta=df_meta.append({'fusion_junc':circ_id,'length':circ_exon_len,'gene_symbol':gene_symbol,
                                'num_miRNA_hits':len(miR_list),'miRNA_list':', '.join(miR_list),
                                'max_seed_matches':max_match,'max_miRNA_list':max_miR_IDs},ignore_index=True)
        
    #add 5 more columns: split reads, chimeric reads, chom, start and end coordinates
    df_meta_con=df_meta.merge(junc_symbol, on=['fusion_junc','gene_symbol'])
    
    df_meta_con.to_csv(out_file)


               