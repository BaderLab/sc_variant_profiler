import subprocess
import re
import os
import glob
import sys
import multiprocessing
import argparse
from Process import Process
from VariantCaller import VariantCaller



def main():
    parser = argparse.ArgumentParser(prog="SC_Mutation" , description="Analyze Mutational properites of  Single Cell RNA-seq data"
                                     , epilog='Enjoy the program! :)')
    parser.add_argument("-A", "--annovar", help="run annovar on my project", dest='function_ann', action='store_true')
    parser.add_argument("-cb", "--cbsniffer", help="Run CbSniffer on the genes", dest='function_cbs', action='store_true')
    parser.add_argument("-eg", "--extract_genes", help="extract the giver gene list", dest='function_eg', action='store_true')
    parser.add_argument("project",  help="the project to run these on ",  type=str)
    parser.add_argument("-gi", help="the gene list input to extract ", type = str , action='store', required=False )
    parser.add_argument("-vc",  help="run variant calling, please choose your tool : Samtools, GATK ", type = str , action='store', required=False )
    parser.add_argument("-sb", "--split-bam", help="Split bamfile with cell barcodes", dest='spliter', action='store_true')
    args = parser.parse_args()
    gene_list = []
    if args.gi:
        gene_list = [line.rstrip('\n') for line in open("genelist")]
    project = args.project
    root = "/home/aiden/gary/Projects/"
    log = open("pipeline.txt", "w+")
    proc = Process(project, gene_list,log, annovar = args.function_ann, cbsniffer = args.function_cbs,
                   extract_genes = args.function_eg)
    if args.vc:
        caller = args.vc
        call = VariantCaller(project, caller ,log)
        if caller == "Samtools":
            call.Samtools()
        elif caller == "GATK":
            call.GATK()


    if args.spliter:
        proc.split_Cells()



    proc.Annovar()
    proc.Extract_genes()
    proc.Cb_sniffer()
if __name__ == '__main__':
    main()





