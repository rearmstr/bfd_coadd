#! /usr/bin/env python

import nsim
import galsim
import bfd
import numpy as np
import coaddsim
import argparse
import multiprocessing
import time
from bfd_coadd import BfdObs


def worker(weight_n,weight_sigma,sigma_step,sigma_max,xy_max,sn_min,ngal,target,template,file,name,label,
           factor,output_dir,index,nobs_cov,use_noise_ps,sigma,psf_seed):

    wname = multiprocessing.current_process().name
    seed=int(np.random.rand()*100000000)+index
    weight = bfd.KSigmaWeight(weight_n, weight_sigma)

    table_multi = None
    table_coadd = None

    sims = nsim.sime.Sim('%s'%(file), seed)

    # Set noise higher and shear to zero and recreate sims
    if template:
        new_config={}
        for key,value in sims.items():
            new_config[key] = value

        if args.sigma is not None:
            new_config['images']['noise'] = sigma/factor
        else:
            new_config['images']['noise'] = sims['images']['noise']/factor

        new_config['shear'] = [0.0, 0.0]
        sims = nsim.sime.Sim(new_config,seed)

    obs_test = sims(psf_seed=psf_seed)
    bfd_test = BfdObs(obs_test, weight, id=0, nda=1./ngal)
    cov_test = bfd_test.moment.get_covariance()
    sigma_flux = np.sqrt(cov_test[0][0,0])
    sigma_xy = np.sqrt(cov_test[1][0,0])

    # run over a number of different coadds to average over different shifts
    sigma_fluxes = []
    sigma_xys = []
    for i in range(nobs_cov):
        obs_test = sims(psf_seed=psf_seed)
        coadd_image = coaddsim.CoaddImages(obs_test, interp='lanczos3')
        coadd = coadd_image.get_mean_coadd()
        bfd_coadd = BfdObs(coadd, weight, id=0, nda=1./ngal, compute_noise_ps=use_noise_ps)

        cov_coadd_test = bfd_coadd.moment.get_covariance()

        sigma_fluxes.append(np.sqrt(cov_coadd_test[0][0,0]))
        sigma_xys.append(np.sqrt(cov_coadd_test[1][0,0]))

    sigma_coadd_flux = np.mean(sigma_fluxes)
    sigma_coadd_xy = np.mean(sigma_xys)


    if template:
        # use errors from noisy measurements
        sigma_flux *= factor
        sigma_xy *= factor
        table_multi = bfd.TemplateTable(weight_n, weight_sigma, sn_min, sigma_xy, sigma_flux,
                                        sigma_step, sigma_max)
        sigma_coadd_flux *= args.factor
        sigma_coadd_xy *= args.factor
        table_coadd = bfd.TemplateTable(weight_n, weight_sigma, sn_min, sigma_coadd_xy, sigma_coadd_flux,
                                    sigma_step, sigma_max)

    else:
        table_multi = bfd.TargetTable(weight_n, weight_sigma, cov=None)
        table_coadd = bfd.TargetTable(weight_n, weight_sigma, cov=None)


    for i in range(ngal):
        if i%(ngal/10)==0 and i>0:
            print wname,"%d%% done"% int(100.0*i/ngal)
        obs_list = sims(psf_seed=psf_seed)
        coadd_image = coaddsim.CoaddImages(obs_list, interp='lanczos3')
        coadd = coadd_image.get_mean_coadd()

        bfd_multi = BfdObs(obs_list, weight, id=i, nda=1./ngal)
        bfd_coadd = BfdObs(coadd, weight, id=i, nda=1./ngal, compute_noise_ps=args.use_noise_ps)

        if template:
            templates = bfd_multi.moment.make_templates(sigma_xy, sn_min=sn_min,
                                                        sigma_flux=sigma_flux, sigma_step=sigma_step,
                                                        sigma_max=sigma_max, xy_max=xy_max)
            for tmpl in templates:
                if tmpl is not None:
                    table_multi.add(tmpl)

            templates = bfd_coadd.moment.make_templates(sigma_coadd_xy, sn_min=sn_min,
                                                        sigma_flux=sigma_coadd_flux, sigma_step=sigma_step,
                                                        sigma_max=sigma_max, xy_max=xy_max)
            for tmpl in templates:
                if tmpl is not None:
                    table_coadd.add(tmpl)
        else:
            xyshift, error, msg = bfd_multi.moment.recenter()
            if error:
                table_multi.addLost()
            else:
                table_multi.add(bfd_multi.moment.get_moment(0,0), xy=xyshift, id=i,
                                covgal=bfd_multi.moment.get_covariance())

            xyshift, error, msg = bfd_coadd.moment.recenter()
            if error:
                table_coadd.addLost()
            else:
                table_coadd.add(bfd_coadd.moment.get_moment(0,0), xy=xyshift, id=i,
                                covgal=bfd_coadd.moment.get_covariance())

    if template:
        table_multi.save('%s/%s_template_multi%s.fits'%(output_dir,name,label))
        table_coadd.save('%s/%s_template_coadd%s.fits'%(output_dir,name,label))
    else:
        table_multi.save('%s/%s_multi%s.fits'%(output_dir,name,label))
        table_coadd.save('%s/%s_coadd%s.fits'%(output_dir,name,label))




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run single file')
    parser.add_argument('--weight_n',default=4,type=int, help='nfor weight function')
    parser.add_argument('--weight_sigma',default=1.,type=np.float64, help='sigma for weight function')
    parser.add_argument('--sigma_step',default=1.2,type=float, help='sigma step for prior')
    parser.add_argument('--sigma_max',default=5.5,type=float, help='sigma max')
    parser.add_argument('--xy_max',default=4,type=float, help='xy max')
    parser.add_argument('--sn_min',default=5,type=float, help='min s/n')
    parser.add_argument('--ngal',default=1000,type=int, help='sigma max')
    parser.add_argument('--target', help='Make target galaxies -OR-', action='store_const', const=True)
    parser.add_argument('--template', help='Make template galaxies', action='store_const', const=True)
    parser.add_argument('--file',default='bdk',type=str, help='sigma max')
    parser.add_argument('--name',default='bdk',type=str, help='sigma max')
    parser.add_argument('--factor',default=100,type=float, help='noise factor for template')
    parser.add_argument('--sigma',default=None,type=float, help='Use this sigma instead of value from file')
    parser.add_argument('--output_dir',default='.', help='output directory')
    parser.add_argument('--start',default=0,type=int, help='output directory')
    parser.add_argument('--njobs',default=48,type=int, help='output directory')
    parser.add_argument('--nobs_cov',default=20, type=int,
                            help='number of observations to compute average covariance on coadd')
    parser.add_argument('--use_noise_ps', dest='use_noise_ps', default=False, action='store_true')
    parser.add_argument('--psf_seed',default=-1,type=int, help='use this seed')

    args = parser.parse_args()
    print args.start

    jobs = []
    for i in range(args.njobs):
        print 'starting',i
        label = '_%d'%(i+args.start)
        arg=(args.weight_n,args.weight_sigma,args.sigma_step,args.sigma_max,args.xy_max,args.sn_min,args.ngal,
             args.target,args.template,args.file,args.name,label,args.factor,args.output_dir,i+args.start,
             args.nobs_cov,args.use_noise_ps,args.sigma,args.psf_seed)
        p = multiprocessing.Process(target=worker, args=arg)
        jobs.append(p)

    for j in jobs:
        j.start()

    for j in jobs:
        j.join()
