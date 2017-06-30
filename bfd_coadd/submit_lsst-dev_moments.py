#! /usr/bin/env python

import time,os,argparse
from collections import deque
import subprocess as sub
import os
import re

parser = argparse.ArgumentParser(description='Run single file')
parser.add_argument('--max_jobs',default=300,type=int,
                    help='number of jobs to run concurrently')
parser.add_argument('--njobs', default=50,type=int,help='number of jobs total to run')
parser.add_argument('--exe',default='python generate_moments.py',
                    help='executable')
parser.add_argument('--hours',type=int,default=3,
                    help='how many hours')
parser.add_argument('--mins',type=int,default=25,
                    help='how many minutes')
parser.add_argument('--arg',default='',
                    help='generic arguments')
parser.add_argument('--ngal', default=200000,type=int,help='dependency')
parser.add_argument('--file',help='file for config')
parser.add_argument('--name',help='name for output')
parser.add_argument('--template', help='Make template galaxies', action='store_const', const=True)
parser.add_argument('--output_dir', default='.',help='Output_directory')
parser.add_argument('--start',default=0, type=int,
                    help='start labels with this number')
parser.add_argument('--jpn',default=48, type=int,
                    help='jobs per node')

args = parser.parse_args()

if os.path.exists(args.output_dir) is False:
    os.makedirs(args.output_dir)
    os.makedirs(args.output_dir+"/submit_files")
# convert options to dictionary to use with format
dict=vars(args)


for ijob in range(args.njobs):

    while True:
        pipe = sub.Popen(['squeue','-u','rea3'],stdout=sub.PIPE)
        # count the number of jobs currently running
        q_out=pipe.communicate()[0]
        num=len(q_out.split('\n'))-1
        if(num<args.max_jobs):break
        time.sleep(0.25)

    name='%s.%d'%(args.name,ijob)
    dict['fullname']=name
    submit_text="""#!/bin/bash
#SBATCH -N 1
#SBATCH -c 48
#SBATCH -J {fullname}
#SBATCH -t {hours}:{mins:02d}:00\n""".format(**dict)

    ilabel = ijob*args.jpn

    if args.template:
            output='%s/log.%s.template.%d'%(args.output_dir,args.name,ilabel)
    else:
            output='%s/log.%s.%d'%(args.output_dir,args.name,ilabel)

    dict['output']=output
    submit_text += 'time python -u generate_mp_moments.py  --ngal %d --start %s --file %s --output_dir %s --name %s --njobs %d' %(args.ngal,ilabel,args.file, args.output_dir,args.name,args.jpn)
    if args.template:
            submit_text += ' --template'
    #submit_text += '>& %s\n'%output

    if args.template:
        submit_file='%s/submit_files/submit_%s_template_%d.cmd' % (args.output_dir,args.name,ijob)
    else:
        submit_file='%s/submit_files/submit_%s_%d.cmd' % (args.output_dir,args.name,ijob)
    ofile=open(submit_file,'w')
    ofile.write(submit_text)
    ofile.close()
    os.system('sbatch %s' % (submit_file))
