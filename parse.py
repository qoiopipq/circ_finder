import sys
import re
from pybedtools import BedTool
import io


#sys inputs:
ref_gene=sys.argv[1]#input reference file: mm10refGene.gtf
junc_file=sys.argv[2]#back_spliced_junction.bed
chim_bam=sys.argv[3]#chimeric.bam


#sys output:
junc_out=sys.argv[4]#parsed_junction.bed
chim_counts=sys.argv[5]#chimer_counts.bed 


def parse(ref_gene, junc_file, chim_bam, junc_out, chim_counts):
    with open(junc_out,'w') as out:
        anno=BedTool(ref_gene).filter(lambda x:x[2]=='exon')
        junc_100k=BedTool(junc_file).filter(lambda f: len(f)<=100000)#junc_100k.bed based on unsorted BED
        intersect_out=BedTool(junc_100k).intersect(anno, wa=True, wb=True, bed=True)#100k_refseq_test.bed
    
        #Extract unique annotated junctions
        junc_anno_uniq=[] 
        for line in intersect_out:
                line=re.split('\t|;|"',str(line))
                #select relevant columns
                if [line[i] for i in [0,1,2,3,4,5,15,12]]not in junc_anno_uniq:
                    junc_anno_uniq.append([line[i] for i in [0,1,2,3,4,5,15,12]])
                    #print("\t".join(map(str,(line[i] for i in [0,1,2,3,4,5,15,12])))+'\n', junc_parsed) 
                    #junc_parsed+=("\t".join(map(str,(line[i] for i in [0,1,2,3,4,5,15,12])))+'\n')
                    print("\t".join(map(str,(line[i] for i in [0,1,2,3,4,5,15,12]))), file=out)
    
    #Bedtools converage
    chim=BedTool(chim_bam)#chimeric.bam
    junc_100k=BedTool(junc_file).filter(lambda f: len(f)<=100000)#junc_100k.bed
    BedTool(junc_100k).coverage(chim).saveas(chim_counts)#chimeric_counts
    

