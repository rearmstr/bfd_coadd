#! /usr/bin/env python

import time,os,argparse
from collections import deque
import subprocess as sub
import re

parser = argparse.ArgumentParser(description='Run single file')
parser.add_argument('--max_jobs',default=300,type=int,
                    help='number of jobs to run concurrently')
parser.add_argument('--njobs', default=50,type=int,help='number of jobs total to run')
parser.add_argument('--name',default='name',
                    help='name of submit job')
parser.add_argument('--exe',default='generate_mp_moments.py',
                    help='executable')
parser.add_argument('--hours',type=int,default=3,
                    help='how many hours')
parser.add_argument('--mins',type=int,default=25,
                    help='how many minutes')
parser.add_argument('--arg',default='',
                    help='generic arguments')
parser.add_argument('--ngal', default=200000,type=int,help='dependency')
parser.add_argument('--file', default='target',help='dependency')
parser.add_argument('--type', default='target',help='dependency')
parser.add_argument('--start',default=0, type=int,
                    help='start labels with this number')
parser.add_argument('--jpn',default=48, type=int,
                    help='jobs per node')
parser.add_argument('--output_dir', default='.',help='Output_directory')
parser.add_argument('--flat_wcs', dest='flat_wcs', default=False, action='store_true')
parser.add_argument('--wait',default=0.25,type=float,
                    help='time between checks')
parser.add_argument('--noise_n',default=48,type=int, help='size of noise image')
args = parser.parse_args()

# convert options to dictionary to use with format
dict=vars(args)

for ii in range(args.njobs):
    ilabel = ii*args.jpn
    while True:
        pipe = sub.Popen(['squeue','-u','rea3'],stdout=sub.PIPE)
        # count the number of jobs currently running
        q_out=pipe.communicate()[0]
        num=len(str(q_out).split('\\n'))-1
        if(num<args.max_jobs):break
        time.sleep(args.wait)
    output='log.%s.%s.%d'%(args.file,args.type, ilabel)
    dict['output']=output
    use_arg = dict['arg']

    use_arg += ' --ngal %d --start %d --file %s --output_dir %s --name %s --njobs %d  --noise_n %d'%(args.ngal,ilabel,args.file,args.output_dir,args.name,args.jpn,args.noise_n)
    #use_arg += ' --ngal %d --start %s --file %s --type %s --njobs 48' %(args.ngal,ilabel,args.file,args.type)
    if args.flat_wcs:
            submit_text += ' --flat_wcs'

    dict['use_arg'] = use_arg
    submit_text="""#!/bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --output={output}
#SBATCH -t {hours}:{mins:02d}:00
{exe} {use_arg} """.format(**dict)


    submit_file='submit_files/submit_%s_%s_%d.cmd' % (args.file,args.type,ilabel)
    ofile=open(submit_file,'w')
    ofile.write(submit_text)
    ofile.close()
    #os.system('sbatch %s' % (submit_file))
