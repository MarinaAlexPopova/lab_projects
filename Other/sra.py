# !/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.04.2022
# @author: Marina Popova
# @contact: marinaalexpopova@yandex.ru

import subprocess

# samples correspond to Het_1, Het_2, Imm_1, Imm_2
sra_numbers = ["SRR18796004", "SRR18796007", "SRR18796011", "SRR18795866"]

# this will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
#for sra_id in sra_numbers:
#    print("Currently downloading: " + sra_id)
#    prefetch = "prefetch " + sra_id
#    print("The command used was: " + prefetch)
#    subprocess.call(prefetch, shell=True)

# this will extract the .sra files from above into a folder named 'fastq'
for sra_id in sra_numbers:
    print("Generating fastq for: " + sra_id)
    fastq_dump = (
        "fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip ./"
        + sra_id + "/" 
         + sra_id 
         + "./"
    )
    print("The command used was: " + fastq_dump)
    subprocess.call(fastq_dump, shell=True)
    
#"fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 " + sra_id
