import os
import glob
import subprocess
from multiprocessing import Pool
from subprocess import Popen, PIPE


class VariantCaller():
    def __init__(self,run_project, Caller ,log):
        self.root = "/home/aiden/hdd/backup_27_04_2020/gary/projects/GBM/unsorted/g6/bam/"
        self.reference_genome = "~/gary/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa"
        self.project = run_project
        self.caller = Caller
        self.log = log

    def GATK(self):
        gref = self.reference_genome
        project_path = self.root + self.project
        print("the project path is :", project_path)
        os.chdir(project_path)
        for file in glob.iglob('**/*.bam', recursive=True):
            print("Kiir", file)
            bam_inp = file[:-4]
            print(bam_inp)
            p1 = Popen(["samtools", "sort", "%s.bam" % (bam_inp),"-o", "%s_samsorted.bam" % (bam_inp)], stdout=PIPE)
            p1.communicate()
            p1.wait()
            p2 = Popen(["gatk", "AddOrReplaceReadGroups", "--INPUT",
                 "%s_samsorted.bam" % (bam_inp), "--OUTPUT", "%s_sort_addgroup.bam" % (bam_inp), "--RGLB", "lib1", "--RGPL","illumina",
                 "--RGPU", "unit1", "--RGSM", "20", "--SORT_ORDER", "coordinate"], stdout=PIPE)
            p2.communicate()
            p2.wait()
            p3 = Popen(["gatk", "MarkDuplicates", "--INPUT","%s_sort_addgroup.bam" % (bam_inp), "--OUTPUT",
                        "%s_sort_addgroup_rawdedupped.bam" % (bam_inp),"--METRICS_FILE", "%s_out.metrics" % (bam_inp)], stdout=PIPE)
            p3.communicate()
            p3.wait()
            bamin = bam_inp+"_sort_addgroup_rawdedupped.bam"
            bamout = bam_inp+"sort_addgroup_dedupped_splitN.bam"
            cmd = "gatk SplitNCigarReads -R ~/gary/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa -I %s -O %s " %(bamin,bamout)
            subprocess.call(cmd, shell=True)
            vcfout = bam_inp+ "_gatk_wholes_genome.vcf"
            cmdl = " gatk HaplotypeCaller -R %s -I %s --dont-use-soft-clipped-bases  --standard-min-confidence-threshold-for-calling 20.0 -O %s" %(gref,bamout,vcfout)
            subprocess.call(cmdl, shell=True)

    def Samtools(self):
        project_path = self.root + self.project
        print("the project path is :", project_path)
        os.chdir(project_path)
        log = self.log
        for file in glob.iglob('**/*.bam', recursive=True):
            inp = file
            outp = file + "_Samtools.vcf"
            print(inp)
            log.write("running samtools on the bam files\n")
            cmd = "bcftools mpileup -Q 30 -A -x -Ou -f ~/gary/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa %s | bcftools call -mv >  %s " % (inp , outp)
            subprocess.call(cmd, shell=True)
            log.write("%s\n" % (cmd))







