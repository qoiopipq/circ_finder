import pandas as pd
import sys

#input files:
converted=sys.argv[1]#UCSC LiftOver converted bed 
target=sys.argv[2]#candidate junctions bed
match=sys.argv[3]#the output file from filter_junc step 

#out file:
out_file=sys.argv[4]#conserved junctions between species 

def conserve(converted, target, match, out_file):
    #set categories
    categories=['both sites', 'homologous','5 prime match','3 prime match', 'no homologous']

    #unique list of mouse circ
    mmu_circ=[]

    #set a data frame 
    df=pd.DataFrame({'chrom':[],'start':[],'end':[],'human_circ':[],'human_split_reads':[],'mouse_circ':[],'mouse_split_reads':[],'categories':[]})

    #mouse circRNAs that converted to human hg38 coordinates
    with open(converted) as mm_hg:
        #human circRNA candidates
         for line in mm_hg:
            line=line.rstrip().split()
            #accommondate +/-2bp window size
            starts=list(range(int(line[1])-2, int(line[1])+3))
            ends=list(range(int(line[2])-2, int(line[2])+3))
  
            with open(target) as hg_circ:
                for circ in hg_circ:
                    circ=circ.rstrip()
                    c=circ.split()
              
                    if line[0]==c[0]:
                        #Perfect match:
                        if any([x in starts[2:3] for x in [int(i) for i in c[1:2]]]) & any([y in ends[2:3] for y in [int(i) for i in c[2:3]]]):              
                            #print(circ, 'mmu: '+line[3], categories[0])
                            mmu_circ.append(line[3])
                            df=df.append({'chrom':c[0], 'start':c[1],'end':c[2],
                                      'human_circ':c[3].split('/')[0], 'human_split_reads':c[3].split('/')[1],
                                      'mouse_circ': line[3].split('/')[0],'mouse_split_reads':line[3].split('/')[1],
                                      'categories':categories[0]},ignore_index=True)
                    
                        #Homologous within +/- 2 nt:
                        if any([x in [starts[i] for i in [0,1,3,4]] for x in [int(i) for i in c[1:2]]]) & any([y in [ends[i] for i in [0,1,3,4]] for y in [int(i) for i in c[2:3]]]):              
                            #print(circ, 'mmu: '+line[3],categories[1])
                            mmu_circ.append(line[3])
                            df=df.append({'chrom':c[0], 'start':c[1],'end':c[2],
                                      'human_circ':c[3].split('/')[0], 'human_split_reads':c[3].split('/')[1],
                                      'mouse_circ': line[3].split('/')[0],'mouse_split_reads':line[3].split('/')[1],
                                      'categories':categories[1]},ignore_index=True)
                    
                    
                        #5' site:    
                        if any([x in starts[2:3] for x in [int(i) for i in c[1:2]]]) & any([y for y in ends[2:3] if y not in [int(i) for i in c[2:3]]]):              
                            #print(circ, 'mmu: '+line[3], categories[2])
                            mmu_circ.append(line[3])
                            df=df.append({'chrom':c[0], 'start':c[1],'end':c[2],
                                      'human_circ':c[3].split('/')[0], 'human_split_reads':c[3].split('/')[1],
                                      'mouse_circ': line[3].split('/')[0],'mouse_split_reads':line[3].split('/')[1],
                                      'categories':categories[2]},ignore_index=True)
                    
                        #3' site:    
                        if any([x for x in starts[2:3] if x not in [int(i) for i in c[1:2]]]) & any([y in ends[2:3] for y in [int(i) for i in c[2:3]]]):              
                            #print(circ, 'mmu: '+line[3], categories[3])
                            mmu_circ.append(line[3])
                            df=df.append({'chrom':c[0], 'start':c[1],'end':c[2],
                                      'human_circ':c[3].split('/')[0], 'human_split_reads':c[3].split('/')[1],
                                      'mouse_circ': line[3].split('/')[0],'mouse_split_reads':line[3].split('/')[1],
                                      'categories':categories[3]},ignore_index=True)
                    
                        
                #no homologous:
                if line[3] not in mmu_circ:
                    df=df.append({'chrom':' ', 'start':' ','end':' ',
                                      'human_circ':' ', 'human_split_reads':' ',
                                      'mouse_circ': line[3].split('/')[0],'mouse_split_reads':line[3].split('/')[1],
                                      'categories':categories[4]},ignore_index=True)
    
    miRNA_match=pd.read_csv(match,header=0,dtype='unicode')#dtype='unicode' to solve the warning msg
    miRNA_match=miRNA_match.loc[:,~miRNA_match.columns.str.contains('^Unnamed')]
    con_match=pd.merge(miRNA_match, df, left_on='fusion_junc', right_on='human_circ')
    
    con_match.to_csv(out_file)
            