import sys
import re
from pybedtools import BedTool
import io

#inputs:
junc_parsed=sys.argv[1]#parsed junctions from parse
known_circ=sys.argv[2]#known_circ.txt from CIRCexplorer2
chim_counts=sys.argv[3]#chimeric counts from parse 
#output:
out_file=sys.argv[4]#circRNAs_reads_symbols.txt


def reannotate(junc_parsed,known_circ,chim_counts, out_file):
    #parsing circRNAs
    #known_circ.txt is the output file (circularRNA_known.txt)from CIRCexplorer2 annotate:
    with open(known_circ) as known_circRNA:
        known_circRNA_loc=[]
        for junc in known_circRNA:
            if re.search("^chr(M|Un_|X_|Y_|\d*_)",junc):continue
            junc = junc.split()
            #start, end, and gene name, according to circularRNA_known.txt format
            known_circRNA_loc.append([*map(str,[junc[i] for i in [1,2,14]])])
            #print ("\t".join(map(str,junc_info)))
 
    circ_as=[]
    categories=["Novel",">10bp"]
    #with open(junc_parsed) as junc_anno_uniq_file:
    with open(junc_parsed) as junc_parsed:
        for line in junc_parsed:
            if re.search("^chr(M|Un_|X_|Y_|\d*_)",line):continue
            line=line.rstrip().split()
            if int(line[2])-int(line[1])<200: continue #exclude back-splicing junctions <200bp
            #need to accommondate +-10bp window size 
            starts=list(range(int(line[1])-10, int(line[1])+11))
            ends=list(range(int(line[2])-10, int(line[2])+11))
       
            for lst in known_circRNA_loc:
                #only check junction start for annotation start, junction end for annotation end 🌟
                #if any(x in starts for x in [int(lst[0])])& any(y in ends for y in [int(lst[1])]) & (line[6] in lst[2]):‼️
                if any([x in starts for x in [int(i) for i in lst[:2]]])& any([y in ends for y in [int(i) for i in lst[:2]]]) & (line[6] in lst[2]):
                    if line+lst[:2] not in circ_as:#[line[i] for i in [0,1,2,3,4,5,15,12]]
                        circ_as.append(line+lst[:2])
                        #print ("\t".join(map(str,line+lst[:2])))
                        break
                
            #Novel circRNAs
            if line+lst[:2] not in circ_as:
                circ_as.append(line+categories[:2])
                #print ("\t".join(line+categories[:2]))

        #Retrieve chimeric reads 
        filtered_chimeric_junc=[]
        #put fusionjunc and chimeric read counts into a dict
        fusion_reads=dict()
        with open(chim_counts) as chim_counts:
            for line in chim_counts:
                #exclude those start with chrM or chrUn_, chrX_, chrY_ 
                if re.search("^chr(M|Un_|X_|Y_|\d*_)",str(line)):continue
                line=re.split('\t|/',str(line))
                #update chimeric read counts: col7-# of back-splicing reads
                line[7]=int(line[7])-int(line[4])
                #exclude those >100k back-splicing junctions
                if int(line[2])-int(line[1])>100000:continue
                line[3]=line[3]+'/'+line[4]
                filtered_chimeric_junc.append(line[:8])
                #print ("\t".join(map(str,line[:8])))
                fusion_reads[line[3]]=line[7]

    
    with open(out_file,'w') as circ_symbols:
        #list for circRNAs on both strands with chimeric read counts and host gene symbols
        for row in circ_as:
            for k,v in fusion_reads.items():
                if k!=row[3]: continue
                row[4]=v
                print("\t".join(map(str,row)), file=circ_symbols)         