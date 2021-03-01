import os
import glob
import subprocess
# import pysam


class Process():
    def __init__(self,run_project,genelist,log , annovar, cbsniffer, extract_genes ):
        self.root ="/home/aiden/hdd/backup_27_04_2020/gary/projects/GBM/unsorted/"
        self.project = run_project
        self.genelist = genelist
        self.cbsniffer = cbsniffer
        self.annovar = annovar
        self.gene_extraction = extract_genes
        self.log = log
        self.files = {}

    def Annovar(self):
        if self.annovar:
            project_path = self.root + self.project
            print("the project path is :", project_path)
            os.chdir(project_path)
            files = []
            for file in glob.iglob('**/*.vcf', recursive=True):
               fname = file.split("/")[-1]
               #get the avinpu
               avinputo = file+".avinput"
               log = self.log
               log.write("making the annovar inputs\n")
               cmd = "convert2annovar.pl -format vcf4 %s > %s " %(file, avinputo)
               subprocess.call(cmd, shell=True)
               log = self.log
               log.write("%s\n" %(cmd))
               annout = file+".annout"
               cmd = "annotate_variation.pl -out %s -build hg38 %s  ~/softwares/annovar/humandb " % (annout, avinputo)
               subprocess.call(cmd, shell=True)
               log.write("%s\n" % (cmd))


    def extractit(self, infile,file):
        genelist = self.genelist
        lfname = file
        def outputer(inp, filename):
            fn = lfname + "." +filename + "_Var.vcf"
            with open(fn, "w") as tst:
                line = "chrm\tstart\tstop\tref\tvar\tgene_name\ttrv_type\n"
                tst.write(line)
                for l in inp.values:
                    inf = l[3:8]
                    inf = "\t".join(str(x) for x in inf)
                    eff = inf + "\t" + l[11] + "\t" + l[1].split()[0]+"\n"
                    tst.write(eff)
            log = self.log
            log.write("the genes are extracted and saved in _Var.vcf files \n")
        for g in genelist:
            filename = g
            s = infile[(infile.gene == g) & ((infile.effect == "nonsynonymous SNV") | (infile.effect == "frameshift insertion") | (
                        infile.effect == "frameshift deletion") | (infile.effect == "frameshift block substitution") | (infile.effect == "stopgain") |
                                             ( infile.effect == "stoploss"))]
            outputer(s, filename)



    def Extract_genes(self):
        if self.gene_extraction:
            project_path = self.root + self.project
            print("the project path is :" , project_path)
            os.chdir(project_path)
            files = []
            for file in glob.iglob('**/*.exonic_variant_function', recursive=True):
                print("the file to process is:", file)
                infile = pd.read_csv(file ,sep="\t", header=None)
                infile.columns = ["redundant", "effect", "gene_info" , "chrm", "start", "end", "ref", "alt", "genotype", "red1", "red2"]
                infile[["gene", "info"]] = infile.gene_info.str.split(":",n = 1, expand= True,)
                self.extractit(infile,file)


    def Cb_sniffer(self):
        if self.cbsniffer:
            project_path = self.root + self.project
            print("the project path is :", project_path)
            os.chdir(project_path)
            for bfile in glob.iglob('**/*codes.tsv.gz', recursive=True):
                cmd = "gunzip %s" %bfile
                subprocess.call(cmd, shell=True)
                print("kard")
            for file in glob.iglob('**/*_.vcf', recursive=True):
                varfile = file
                path = varfile.split("/")
                path2 = "/".join(e for e in path[0:2])
                barcodes = path2 + "/barcodes/barcodes.tsv"
                bam = path2 + "/*.bam"
                print(bam)
                cmd = "python3 cb_sniffer.py %s %s %s %s"  %(bam,varfile,barcodes, varfile)
                subprocess.call(cmd, shell=True)
                print("don't worry its doing its job :D")

    def split_Cells(self):
        project_path = self.root + self.project
        print("the project path is :", project_path)
        os.chdir(project_path)
        for bsfile in glob.iglob('**/*.bamsorted.bam', recursive=True):
            ### Written by https://github.com/herrinca
            print(bsfile)
            unsplit_file = bsfile
            out_dir = "".join(bsfile.split("/")[:-1]) + "/Cells_bamFiles/"
            os.mkdir(out_dir)
            # variable to hold barcode index
            CB_hold = 'unset'
            itr = 0
            # read in upsplit file and loop reads by line
            samfile = pysam.AlignmentFile(unsplit_file, "rb")
            for read in samfile.fetch(until_eof=True):
                # barcode itr for current read
                CB_itr = read.get_tag('BC')
                # if change in barcode or first line; open new file
                if (CB_itr != CB_hold or itr == 0):
                    # close previous split file, only if not first read in file
                    if (itr != 0):
                        split_file.close()
                    CB_hold = CB_itr
                    itr += 1
                    split_file = pysam.AlignmentFile(out_dir + "CB_{}.bam".format(itr), "wb", template=samfile)

                # write read with same barcode to file
                split_file.write(read)
            split_file.close()
            samfile.close()


