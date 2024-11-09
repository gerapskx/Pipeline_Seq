#!/usr/bin/env python3

import os
from glob import glob
import subprocess

# initialization
work_dir = './Cpepo_As'
samples = {
    'Cpepo_100_1_S7_R1_001.fastq': 'SRR31111300',
    'Cpepo_100_2_S8_R1_001.fastq': 'SRR31111299',
'Cpepo_100_3_S9_R1_001.fastq': 'SRR31111298',
'Cpepo_200_1_S10_R1_001.fastq': 'SRR31111297',
'Cpepo_200_2_S11_R1_001.fastq': 'SRR31111296',
'Cpepo_200_3_S12_R1_001.fastq': 'SRR31111295',
'Cpepo_50_1_S4_R1_001.fastq': 'SRR31111305',
    'Cpepo_50_2_S5_R1_001.fastq': 'SRR31111304',
'Cpepo_50_3_S6_R1_001.fastq': 'SRR31111301',
'Cpepo_CTRL_1_S1_R1_001.fastq': 'SRR31111294',
'Cpepo_CTRL_2_S2_R1_001.fastq': 'SRR31111303',
'Cpepo_CTRL_3_S3_R1_001.fastq': 'SRR31111302',
}


# downloading each given file
for sample_id in samples:
    print('Currently downloading: ' + samples[sample_id])

    # downloading/converting the files
    cmd_prefetch = 'prefetch --output-directory {:s} --progress {:s}'.format(work_dir, samples[sample_id])
    print('\trunning: ' + cmd_prefetch)
    subprocess.call(cmd_prefetch, shell=True)

    cmd_fastqdump = 'fastq-dump --outdir {:s} --skip-technical --readids '.format(work_dir) + \
                    '--read-filter pass --dumpbase --split-3 --clip ' + \
                    '{:s}/{:s}/{:s}.sra'.format(work_dir, samples[sample_id], samples[sample_id])
    print('\trunning: ' + cmd_fastqdump)
    subprocess.call(cmd_fastqdump, shell=True)

    # compressing the fastqs
    for fq_name in glob('{:s}/{:s}*.fastq'.format(work_dir, samples[sample_id])):
        cmd_compress = 'gzip -c {:s} > {:s}/{:s}_{:s}.gz'.format(fq_name, work_dir, sample_id, os.path.basename(fq_name))
        print('\trunning: ' + cmd_compress)
        subprocess.call(cmd_compress, shell=True)
        os.remove(fq_name)

    # clean up
    cmd_rmdir = 'rm -r {:s}/{:s}'.format(work_dir, samples[sample_id])
    print('\trunning: ' + cmd_rmdir)
    subprocess.call(cmd_rmdir, shell=True)

