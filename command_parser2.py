import sys
import argparse

#top-level parser
parser=argparse.ArgumentParser(description='circRNA parsing and seed matching.')
subparsers=parser.add_subparsers()

#seed_extract subcommand
parser_extract=subparsers.add_parser('seed_extract',help='extracting seed sequences from miRNAs')
parser_extract.add_argument('mir_seq',  help="miRNA sequences")
parser_extract.add_argument('out_file',  help="seed sequences")

#parse subcommand
parser_parse=subparsers.add_parser('parse',help='parsing back splicing junctions')
parser_parse.add_argument('ref_gene', help="reference file")
parser_parse.add_argument('junc_file',help="back_spliced_junction.bed from CIRCexplorer2")
parser_parse.add_argument('chim_bam', help='chimeric.bam')
parser_parse.add_argument('junc_out',help="parsed_junc.bed")
parser_parse.add_argument('chim_counts', help='chimeric_counts.bed')

#reannotate subcommand
parser_reannotate=subparsers.add_parser('reannotate', help='reannotate circRNAs')
parser_reannotate.add_argument('junc_parsed', help='processed junctions from parse step')
parser_reannotate.add_argument('known_circ',help='knwon circRNAs from CIRCexplorer2')
parser_reannotate.add_argument('chim_counts', help='chimeric read counts')
parser_reannotate.add_argument('out_file', help='circRNAs_reads_symbols.txt')

#filter junctions 
parser_filter=subparsers.add_parser('filter_junc', help='filtering junctions')
parser_filter.add_argument('circ_annotated', help='circRNAs_reads_symbols.txt from reannotate step')
parser_filter.add_argument('junc_file', help='back_spliced_junction.bed from CIRCexplorer2')
parser_filter.add_argument('out_file', help='filtered_junc.bed')
parser_filter.add_argument('out_filtered_df', help='junc_filtered_readcounts.csv')


#extract exonic sequences and seed matching 
parser_match=subparsers.add_parser('match', help='seed matching for selected circRNAs')
parser_match.add_argument('junc_cand', help='filtered_junc.bed from filter_junc step')
parser_match.add_argument('ref_exon', help='reference exons BED')
parser_match.add_argument('genome', help='genome.fa')
parser_match.add_argument('fil_junc', help='filtered_readcounts.csv from filter_junc step')
parser_match.add_argument('seeds', help='miRNA seed sequences')
parser_match.add_argument('out_file', help='circ_seed_matches.csv')
                      

#find conserved junctions between species
parser_conserve=subparsers.add_parser('conserve', help='detecting conserved junctions between species')
parser_conserve.add_argument('converted', help='UCSC LiftOver converted junctions')
parser_conserve.add_argument('target', help='Candidate/Filtered junctions that are predicted for target organism')
parser_conserve.add_argument('match', help='Output file from match step for target organism')
parser_conserve.add_argument('out_file', help='Conserved junctions between species')

args=parser.parse_args()   

def main():
    if len(sys.argv[1:])==0:
        parser.print_help()
    elif sys.argv[1]=='seed_extract':
        from seed_extract import seed_extract
        seed_extract(sys.argv[2], sys.argv[3]) 
    elif sys.argv[1]=='parse': #parsing 
        from parse import parse 
        parse(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1]=='reannotate': #reannotate 
        from reannotate import reannotate
        reannotate(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif sys.argv[1]=='filter_junc': #filter junctions
        from filter_junc import filter_junctions
        filter_junctions(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif sys.argv[1]=='match':# seed matching 
        from match import extract_exon_seq
        extract_exon_seq(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
                          
    elif sys.argv[1]=='conserve': #conserved junctions
        from conserve import conserve
        conserve(sys.argv[2], sys.argv[3], sys.argv[4],sys.argv[5])
   
        
        
        
if __name__ == '__main__':
    main()