#################################################################################################### DEPENDENCIES

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

##################################################################################################### ARGUMENTS

parser = argparse.ArgumentParser(prog='ants lcWGS pipeline', description='ants lcWGS pipeline')

parser.add_argument('-1', '--one',metavar='\b', help='paired reads rx mate')
parser.add_argument('-2', '--two',metavar='\b', help='paired reads lx mate')
parser.add_argument('-u', '--unpaired',metavar='\b', help='single reads')
parser.add_argument('-s', '--subset',default="1000000",metavar='\b', help='subset million of reads for mtgen assembly')
parser.add_argument('-l', '--mtlen',default="1000",metavar='\b', help='mt contigs length cutoff')
parser.add_argument('-c', '--mtcov',default="10",metavar='\b', help='mt contigs cov... cutoff')
parser.add_argument('-d', '--db_uc',metavar='\b', help='ultraconserved elements database')
parser.add_argument('-o', '--output',metavar='\b', help='output folder')
parser.add_argument('-m', '--mtonly',action='store_true', help='output folder')
parser.add_argument('-t', '--threads',metavar='\b', help='number of threads')
parser.add_argument('-v', '--verbose',action='store_false', help='keeps temporary folder and files')
parser.add_argument('-e', '--erase', action='store_true', help='erases and rewrites a pre existing output folder')

args=parser.parse_args()

if args.db_uc == None and args.mtonly == False:
    print("\nWARNING! For UCE reconstruction a bait set must be specified\n")
    quit()

##################################################################################################### FOLDER SYSTEM

print("\n" + datetime.datetime.now().strftime("%H:%M:%S %Y-%m-%d") + " start")

if args.erase == False and path.exists(args.output):
    print("\nWARNING! An output folder with the same name already exists! Use --erase to overwrite \n")
    quit()

if path.exists(args.output):
	shutil.rmtree(args.output + "/tmp")

os.makedirs(args.output + "/tmp")
os.chdir(args.output + "/tmp")

if not args.unpaired:
	mates1 =  "../../" + args.one
	mates2 =  "../../" + args.two
else:
	mates = "../../" + args.unpaired

################################################################################################ TRIM

print(datetime.datetime.now().strftime("%H:%M:%S") + ' reads trimming')

if not args.unpaired:
	with open('log_trimmomatic', 'w') as log_trimmomatic:
		subprocess.call(['trimmomatic', 'PE',
		mates1,mates2,'-baseout','trimmomatic',PE_ADAPTERS_PATH,
		'LEADING:20','TRAILING:20',
		'SLIDINGWINDOW:5:20','MINLEN:95'],stdout=log_trimmomatic,stderr=log_trimmomatic)
else:
	with open('log_trimmomatic', 'w') as log_trimmomatic:
                subprocess.call(['trimmomatic', 'SE',
                mates,'trimmomatic_U',SE_ADAPTERS_PATH,
                'LEADING:20','TRAILING:20',
                'SLIDINGWINDOW:5:20','MINLEN:95'],stdout=log_trimmomatic,stderr=log_trimmomatic)


############################################################################################### QC

print(datetime.datetime.now().strftime("%H:%M:%S") + ' reads QC')

if not args.unpaired:
	with open('log_fastqc', 'w') as log_fastqc:
		subprocess.call(["fastqc", "-o folder_fastqc", "trimmomatic_1P", "trimmomatic_2P"],stdout=log_fastqc,stderr=log_fastqc)
	os.rename('trimmomatic_1P' ,"trimmomatic_1P.fq")
	os.rename('trimmomatic_2P' ,"trimmomatic_2P.fq")
else:
	with open('log_fastqc', 'w') as log_fastqc:
		subprocess.call(["fastqc", "-o folder_fastqc", "trimmomatic_U"],stdout=log_fastqc,stderr=log_fastqc)
	os.rename('trimmomatic_U' ,"trimmomatic_U.fq")

############################################################################################### MT ASSEMBLY

print(datetime.datetime.now().strftime("%H:%M:%S") + ' mt assembly')

if not args.unpaired:
	arg1="in=trimmomatic_1P.fq"
	arg2="in2=trimmomatic_2P.fq"

	if args.subset is not None:
        	arg3="reads=" + args.subset
        	with open('log_reformat', 'w') as log_reformat:
                	subprocess.call(['reformat.sh',arg1,arg2,arg3,'out=1.fq','out2=2.fq'],stdout=log_reformat, stderr=log_reformat)
	if args.subset is None:
        	with open('log_reformat', 'w') as log_reformat:
                	subprocess.call(['reformat.sh',arg1,arg2,'READS=-1','out=1.fq','out2=2.fq'],stdout=log_reformat, stderr=log_reformat)

	with open('log_mt_spades', 'w') as log_mt_spades:
		subprocess.call(["spades.py","-1", "1.fq", "-2", "2.fq",
		"-o","folder_mt_spades", "-t", args.threads, "--cov-cutoff", "off"],stdout=log_mt_spades, stderr=log_mt_spades)

else:
	arg1="in=trimmomatic_U.fq"
	if args.subset is not None:
		arg3="reads=" + args.subset
		with open('log_reformat', 'w') as log_reformat:
			subprocess.call(['reformat.sh',arg1,arg3,'out=U.fq'],stdout=log_reformat, stderr=log_reformat)
	if args.subset is None:
		with open('log_reformat', 'w') as log_reformat:
			subprocess.call(['reformat.sh',arg1,'READS=-1','out=U.fq'],stdout=log_reformat, stderr=log_reformat)

	with open('log_mt_spades', 'w') as log_mt_spades:
		subprocess.call(["spades.py","-s", "U.fq",
		"-o","folder_mt_spades", "-t", args.threads, "--cov-cutoff", "off"],stdout=log_mt_spades, stderr=log_mt_spades)


mt_seqenc = open("mtcontig.fa",'w')

for record in SeqIO.parse("folder_mt_spades/contigs.fasta", "fasta"):
	len = float(record.id.split('_')[3])
	cov = float(record.id.split('_')[5])
	if len > float(args.mtlen) and cov > float(args.mtcov):
		record.id = record.id.split('_length')[0]
		record.description = record.id.split('_length')[0]
		SeqIO.write(record, mt_seqenc, "fasta")
mt_seqenc.close()

############################################################################################### MT ANNOT

print(datetime.datetime.now().strftime("%H:%M:%S") + ' mt annotate')

with open('log_mitoz', 'w') as log_mitoz:
        subprocess.call(["mitoz","annotate","--fastafiles","mtcontig.fa","--outprefix","annot","--thread_number",args.threads,"--species_name",args.output], stdout=log_mitoz, stderr=log_mitoz)

############################################################################################### MT ONLY

if args.mtonly == True:
	os.chdir('..')
	shutil.copy('tmp/annot.mtcontig.fa.result/annot_mtcontig.fa_mitoscaf.fa.gbf','.')
	os.rename('annot_mtcontig.fa_mitoscaf.fa.gbf','mt.gbf')
	shutil.copy('tmp/annot.mtcontig.fa.result/annot_mtcontig.fa_mitoscaf.fa.gbf.fasta','.')
	os.rename('annot_mtcontig.fa_mitoscaf.fa.gbf.fasta','mt.fasta')
	shutil.copy('tmp/annot.mtcontig.fa.result/summary.txt','.')
	os.rename('summary.txt','mt_summary.txt')

	if os.path.getsize('tmp/annot.mtcontig.fa.result/summary.txt') == 0:
		print("\nWARNING! No mt genome assembled! \n")
		quit()

	if args.verbose == False :
        	pass
	else:
        	shutil.rmtree("tmp/", ignore_errors=True)
	print(datetime.datetime.now().strftime("%H:%M:%S") + ' finished - mt only')
	quit()

############################################################################################### NC ASSEMBLY

if args.subset is not None:

	print(datetime.datetime.now().strftime("%H:%M:%S") + ' nc assembly')

	if not args.unpaired:

		with open('log_nc_spades', 'w') as log_nc_spades:
        		subprocess.call(["spades.py","-1", "trimmomatic_1P.fq", "-2", "trimmomatic_2P.fq",
        		"-o","folder_nc_spades", "-t", args.threads, "--cov-cutoff", "off"], stdout=log_nc_spades, stderr=log_nc_spades)

	else:

		with open('log_nc_spades', 'w') as log_nc_spades:
			subprocess.call(["spades.py","-s", "U.fq",
			"-o","folder_nc_spades", "-t", args.threads, "--cov-cutoff", "off"],stdout=log_nc_spades, stderr=log_nc_spades)


############################################################################################## NC ANNOT UCEs

print(datetime.datetime.now().strftime("%H:%M:%S") + ' nc annotate UCEs')

if os.path.getsize('folder_nc_spades/contigs.fasta') == 0:
	print("\nWARNING! No nc genome assembled! \n")
	quit()

uce_db = '../../' + args.db_uc


if args.subset is not None:
	with open('log_phyluce', 'w') as log_phyluce:
        	subprocess.call(["phyluce_assembly_match_contigs_to_probes",
        	"--contigs","folder_nc_spades","--probes",uce_db,"--output","folder_phyluce","--log-path","folder_phyluce"], stdout=log_phyluce, stderr=log_phyluce)
else:
	with open('log_phyluce', 'w') as log_phyluce:
		subprocess.call(["phyluce_assembly_match_contigs_to_probes",
		"--contigs","folder_mt_spades","--probes",uce_db,"--output","folder_phyluce","--log-path","folder_phyluce"], stdout=log_phyluce, stderr=log_phyluce)

############################################################################################# NC ANNOT BUSCO

print(datetime.datetime.now().strftime("%H:%M:%S") + ' nc annotate busco')

if args.subset is not None:
	with open('log_busco', 'w') as log_busco:
		subprocess.call(["busco", "-i", "folder_nc_spades/contigs.fasta","-o", "folder_busco", "--mode", "geno",
		"--cpu", args.threads, "-l", "hymenoptera", "-f"],stdout=log_busco,stderr=log_busco)
else:
	with open('log_busco', 'w') as log_busco:
		subprocess.call(["busco", "-i", "folder_mt_spades/contigs.fasta","-o", "folder_busco", "--mode", "geno", 
		"--cpu", args.threads, "-l", "hymenoptera", "-f"],stdout=log_busco,stderr=log_busco)

############################################################################################ WRAP UP

print(datetime.datetime.now().strftime("%H:%M:%S") + ' wrapping up')

os.chdir('..')

shutil.copy('tmp/annot.mtcontig.fa.result/annot_mtcontig.fa_mitoscaf.fa.gbf','.')
os.rename('annot_mtcontig.fa_mitoscaf.fa.gbf','mt.gbf')
shutil.copy('tmp/annot.mtcontig.fa.result/annot_mtcontig.fa_mitoscaf.fa.gbf.fasta','.')
os.rename('annot_mtcontig.fa_mitoscaf.fa.gbf.fasta','mt.fasta')
shutil.copy('tmp/annot.mtcontig.fa.result/summary.txt','.')
os.rename('summary.txt','mt_summary.txt')

if args.subset is not None:
	shutil.copy('tmp/folder_nc_spades/contigs.fasta','.')
	os.rename('contigs.fasta','nc.fasta')
else:
	os.rename('mt.fasta','assembly.fasta')

shutil.copy('tmp/folder_busco/run_hymenoptera_odb10/short_summary.txt','.')
os.rename('short_summary.txt','nc_summary.txt')
copy_tree('tmp/folder_busco/run_hymenoptera_odb10/busco_sequences/','.')

if args.verbose == False :
	pass
else:
	shutil.rmtree("tmp/", ignore_errors=True)

print(datetime.datetime.now().strftime("%H:%M:%S") + ' finished - mt + nc')
