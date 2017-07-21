#! /usr/bin/env python

import time,os,argparse
from collections import deque
import subprocess as sub
import re
import os

parser = argparse.ArgumentParser(description='Run single file')
parser.add_argument('--max_jobs',default=300,type=int,
                    help='number of jobs to run concurrently')
parser.add_argument('--njobs', default=10,type=int,help='number of jobs total to run')
parser.add_argument('--nfiles', default=100,type=int,help='number of jobs total to run')
parser.add_argument('--exe',default='tableIntegrate',
                    help='executable')
parser.add_argument('--hours',type=int,default=3,
                    help='how many hours')
parser.add_argument('--mins',type=int,default=25,
                    help='how many minutes')
parser.add_argument('--dir', default='.',help='dependency')
parser.add_argument('--label', default='target_target_multi',help='dependency')
parser.add_argument('--template_label', default='target_target_multi',help='dependency')
parser.add_argument('--template', default='template_template_multi.fits',help='dependency')
parser.add_argument('--start',default=0, type=int,
                    help='start labels with this number')
parser.add_argument('--config',default='',
                    help='config file for integration')
parser.add_argument('--sn_bins',default='5,25',
                    help='bins of S/N')
parser.add_argument('--check_exists',default=False,
                    help='check if file already exists')


args = parser.parse_args()
if os.path.exists(args.dir) is False:
    os.makedirs(args.dir+"/submit_files")
# convert options to dictionary to use with format
dict=vars(args)

file_list=[]
for ii in range(args.nfiles):
        ilabel = args.start + ii
        if os.path.exists('%s/%s_%d.fits'%(args.dir,args.label,ilabel)) is False:
            continue
        file_list.append(ilabel)

print 'Found ',len(file_list),' possible files'
file_lists=[file_list[i::args.njobs] for i in range(args.njobs)]

for ijob,file_list in enumerate(file_lists):
    print ijob
    # Wait until the number of jobs is low enough
    while True:
        pipe = sub.Popen(['squeue','-u','rea3'],stdout=sub.PIPE)
        # count the number of jobs currently running
        q_out=pipe.communicate()[0]
        num=len(q_out.split('\n'))-1
        if(num<args.max_jobs):break
        time.sleep(0.25)

    output='%s/log.integrate.%s_%s_%d'%(args.dir, args.label, args.template_label, ijob)
    dict['output']=output

    submit_text="""#!/bin/bash
#SBATCH -N 1
#SBATCH -c 48
#SBATCH --output={output}
#SBATCH -t {hours}:{mins:02d}:00
""".format(**dict)
    print file_list
    for ii in file_list:
        ilabel = args.start + ii

        filename = '%s/pqr_%s_%s_%s.fits'% (args.dir,args.label,args.template_label,ilabel)
        if args.check_exists:
                if os.path.exists(filename):
                    print filename,'already exists'
                    continue

        use_arg = '  %s -targetFile %s/%s_%d.fits -templateFile %s -selectSn=%s -pqrFile %s/pqr_%s_%s_%s.fits -useAveCov True' %(args.config,
                                                                                                                                 args.dir,args.label,ilabel,
                                                                                                                                 args.template,
                                                                                                                                 args.sn_bins,
                                                                                                                                 args.dir,args.label,args.template_label,ilabel)

        dict['use_arg'] = use_arg

        submit_text +="{exe} {use_arg} \n".format(**dict)

    submit_file='%s/submit_files/submit_integrate_%s_%s_%d.cmd' % (args.dir,args.label,args.template_label,ijob)
    ofile=open(submit_file,'w')
    ofile.write(submit_text)
    ofile.close()
    os.system('sbatch %s' % (submit_file))
