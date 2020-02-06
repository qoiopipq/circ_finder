import sys

mir_seq=sys.argv[1]#updated2_hsa_miR.fa
out_file=sys.argv[2]#test_seeds.fa


def seed_extract(mir_seq, out_file):
    with open(out_file,'w') as out:
        with open(mir_seq,'r') as miR:
            for line in miR:
                if line.startswith(">"):
                    h=line.strip().split('>')[1]#header
                    s=list(next(miR).strip())
                    s=['T' if x=='U' else x for x in s]
                    #extract seed sequences: nucleotides 2-7 from mature miRNAs
                    print('>'+h+'\n'+''.join(s[1:7]),file=out)
                

                
  
     