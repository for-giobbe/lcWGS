#################################################################################################### dependencies

import os
import sys
import glob
import shutil
import datetime
import argparse
import subprocess
from os import path
from Bio import SeqIO
from distutils.dir_util import copy_tree

PE_ADAPTERS_PATH="ILLUMINACLIP:/home/PERSONALE/giobbe.forni2/miniconda3/envs/phyluce-1.7.3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:1:30:10"
SE_ADAPTERS_PATH="ILLUMINACLIP:/home/PERSONALE/giobbe.forni2/miniconda3/envs/phyluce-1.7.3/share/trimmomatic-0.39-2/adapters/TruSeq3-SE.fa:1:30:10"

##################################################################################################### arguments

parser = argparse.ArgumentParser(prog='ants lcWGS pipeline', description='ants lcWGS pipeline')

parser.add_argument('-1', '--one', metavar='\b', help='paired reads rx mate')
parser.add_argument('-2', '--two', metavar='\b', help='paired reads lx mate')
parser.add_argument('-u', '--unpaired', metavar='\b', help='single reads')
parser.add_argument('-s', '--subset', metavar='\b', help='subset million of reads')
parser.add_argument('-du', '--database_uc', metavar='\b', help='ultraconserved elements database')
parser.add_argument('-dm', '--database_mt', metavar='\b', help='mitogenomes database')
parser.add_argument('-o', '--output', metavar='\b', help='output folder')
parser.add_argument('-r', '--reference', metavar='\b', help='reference mitochondrial database')
parser.add_argument('-l', '--mtlen', metavar='\b', help='mt contigs length cutoff')
parser.add_argument('-c', '--mtcov', metavar='\b', help='mt contigs cov... cutoff')
parser.add_argument('-t', '--threads', metavar='\b', help='number of threads')
parser.add_argument('-v', '--verbose', action='store_false', help='keeps temporary folder and files')
parser.add_argument('-e', '--erase', action='store_true', help='erases and rewrites a pre existing output folder')

args=parser.parse_args()

##################################################################################################### folder system

print("\n# start on    " + datetime.datetime.now().strftime("%H:%M:%S %Y-%m-%d"))

if args.erase == False and path.exists(args.output):
    print("\n WARNING! An output folder with the same name already exists! Use --erase to overwrite \n")
    quit()

if path.exists(args.output):
	shutil.rmtree(args.output + "/tmp")

os.makedirs(args.output + "/tmp")
os.chdir(args.output + "/tmp")

mates1 =  "../../" + args.one
mates2 =  "../../" + args.two

############################################################################################################ SUBSET

arg1="in=" + mates1
arg2="in2=" + mates2

if args.subset is not None:
	arg3="reads=" + args.subset
	with open('log_reformat', 'w') as reformat_log:
        	subprocess.call(['reformat.sh',arg1,arg2,arg3,'out=1.fq','out2=2.fq'], stdout=reformat_log, stderr=reformat_log)
if args.subset is None:
        with open('log_reformat', 'w') as reformat_log:
                subprocess.call(['reformat.sh',arg1,arg2,'READS=-1','out=1.fq','out2=2.fq'], stdout=reformat_log, stderr=reformat_log)

print("\n# reads trimming    " + datetime.datetime.now().strftime("%H:%M:%S")) #####################################

with open('log_trimmomatic', 'w') as trimmomatic_log:
	subprocess.call(['trimmomatic', 'PE',
	'1.fq', '2.fq', '-baseout', 'trimmomatic', PE_ADAPTERS_PATH,
	'LEADING:20', 'TRAILING:20',
	'SLIDINGWINDOW:5:20', 'MINLEN:95'], stdout=trimmomatic_log, stderr=trimmomatic_log)

print("# reads QC    " + datetime.datetime.now().strftime("%H:%M:%S")) ##############################################

with open('log_fastqc', 'w') as fastqc_log:
	subprocess.call(["fastqc", "-o folder_fastqc", "trimmomatic_1P", "trimmomatic_2P"], stdout=fastqc_log, stderr=fastqc_log)  

os.rename('trimmomatic_1P' ,"trimmomatic_1P.fq")
os.rename('trimmomatic_2P' ,"trimmomatic_2P.fq")

print("# assembly    " + datetime.datetime.now().strftime("%H:%M:%S")) ##############################################

with open('log_spades', 'w') as spades_log:
	subprocess.call(["spades.py","-1", "trimmomatic_1P.fq", "-2", "trimmomatic_2P.fq",
	"-o","folder_spades", "-t", args.threads, "--cov-cutoff", "off"], stdout=spades_log, stderr=spades_log)

mt_seqenc = open("mtcontig.fa",'w')

for record in SeqIO.parse("folder_spades/contigs.fasta", "fasta"):
	len = float(record.id.split('_')[3])
	cov = float(record.id.split('_')[5])
	if len > float(args.mtlen) and cov > float(args.mtcov):
		record.id = record.id.split('_length')[0]
		record.description = record.id.split('_length')[0]
		SeqIO.write(record, mt_seqenc, "fasta")
mt_seqenc.close()

print("#[" + datetime.datetime.now().strftime("%H:%M:%S") + "]	annotate mitogenome") ##############################################

with open('log_mitoz', 'w') as mitoz_log:
        subprocess.call(["mitoz","annotate","--fastafiles","mtcontig.fa","--outprefix","annot","--thread_number",args.threads], stdout=mitoz_log, stderr=mitoz_log)








exit()









print("#    UCE annotation  " + datetime.datetime.now().strftime("%H:%M:%S")) #####################################################################################################################################################

with open('phyluce_busco', 'w') as phyluce_log:
	subprocess.call(["phyluce_assembly_match_contigs_to_probes",'--min_coverage',80,'--min_identity',80,"--contigs","folder_spades/contigs.fasta","--probes",database_uc])

print("#    BUSCO genes annotation  " + datetime.datetime.now().strftime("%H:%M:%S")) #####################################################################################################################################################

with open('log_busco', 'w') as busco_log:
	subprocess.call(["busco", "-i", "folder_spades/contigs.fasta","-o", "folder_busco", "--mode", "geno", "--cpu", args.threads, "-l", "hymenoptera_odb10", "-f"], stdout=busco_log, stderr=busco_log)

print("#    wrapping up       " + datetime.datetime.now().strftime("%H:%M:%S")) ########################################################################################################################$

#os.rename('folder_busco/run_hymenoptera_odb10/short_summary.txt', 'folder_busco/run_hymenoptera_odb10/busco_summary.txt')
#shutil.copy('folder_busco/run_hymenoptera_odb10/busco_summary.txt', '..')
#copy_tree('folder_busco/run_hymenoptera_odb10/busco_sequences/', '..')

#copy_tree('annotate.mitoscaf.hmmtblout.besthit.sim.filtered-by-taxa.fa.result', '..')

os.chdir('..')

if args.verbose == False :
	pass
else:
	shutil.rmtree("tmp/", ignore_errors=True)
