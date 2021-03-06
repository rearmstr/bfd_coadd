#! /usr/bin/env python

import nsim
import galsim
import bfd
import numpy as np
import coaddsim
import argparse
from bfd_coadd import BfdObs


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
parser.add_argument('--noise_range', help='Noise is using range, needed for templates', action='store_const', const=True,default=None)
parser.add_argument('--file',default='bdk',type=str, help='sigma max')
parser.add_argument('--name',default='bdk',type=str, help='sigma max')
parser.add_argument('--label',default='',type=str, help='label')
parser.add_argument('--factor',default=100,type=float, help='noise factor for template')
parser.add_argument('--sigma',default=None,type=float, help='Use this sigma instead of value from file')
parser.add_argument('--output_dir',default='.', help='output directory')
parser.add_argument('--nobs_cov',default=20, type=int,
                    help='number of observations to compute average covariance on coadd')
parser.add_argument('--use_noise_ps', dest='use_noise_ps', default=False, action='store_true')
parser.add_argument('--seed',default=None,type=int, help='use this seed')
parser.add_argument('--psf_seed',default=-1,type=int, help='use this seed')
parser.add_argument('--interp',default='lanczos3',type=str, help='interpolation kernel')
parser.add_argument('--flat_wcs',dest='flat_wcs', default=False, action='store_true')
parser.add_argument('--noise_n',default=48,type=int, help='size of noise image')

args = parser.parse_args()


if args.seed is None:
    seed = np.random.randint(0,2**30)
else:
    seed = args.seed
weight = bfd.KSigmaWeight(args.weight_n, args.weight_sigma)

table_multi = None
table_coadd = None

sims = nsim.sime.Sim('%s'%(args.file), seed)

# Set noise higher and shear to zero and recreate sims
if args.template:
    new_config={}
    for key,value in sims.items():
        new_config[key] = value

    if args.noise_range is None:
        if args.sigma is not None:
            new_config['images']['noise'] = args.sigma/args.factor
        else:
            new_config['images']['noise'] = sims['images']['noise']/args.factor
    else:
        nrange=new_config['images']['noise']['sigma_range']
        new_config['images']['noise']['sigma_range'] = [nrange[0]/args.factor, nrange[1]/args.factor]
    print(new_config)

    new_config['shear'] = [0.0, 0.0]
    sims = nsim.sime.Sim(new_config,seed)

obs_test = sims(psf_seed=args.psf_seed)
bfd_test = BfdObs(obs_test, weight, id=0, nda=1./args.ngal)
cov_test = bfd_test.moment.get_covariance()
sigma_flux = np.sqrt(cov_test[0][0,0])
sigma_xy = np.sqrt(cov_test[1][0,0])

# run over a number of different coadds to average over different shifts
sigma_fluxes = []
sigma_xys = []
sigma_multi_fluxes = []
sigma_multi_xys = []
for i in range(args.nobs_cov):
    obs_test = sims(psf_seed=args.psf_seed)
    bfd_multi = BfdObs(obs_test, weight, id=0, nda=1./args.ngal)
    cov_multi = bfd_multi.moment.get_covariance()
    sigma_multi_flux = np.sqrt(cov_multi[0][0,0])
    sigma_multi_xy = np.sqrt(cov_multi[1][0,0])
    sigma_multi_fluxes.append(sigma_multi_flux)
    sigma_multi_xys.append(sigma_multi_xy)

    coadd_image = coaddsim.CoaddImages(obs_test, interp=args.interp, flat_wcs=args.flat_wcs, noise_n=args.noise_n)
    coadd = coadd_image.get_mean_coadd()
    bfd_coadd = BfdObs(coadd, weight, id=0, nda=1./args.ngal, compute_noise_ps=args.use_noise_ps)

    cov_coadd_test = bfd_coadd.moment.get_covariance()

    sigma_fluxes.append(np.sqrt(cov_coadd_test[0][0,0]))
    sigma_xys.append(np.sqrt(cov_coadd_test[1][0,0]))

sigma_coadd_flux = np.median(sigma_fluxes)
sigma_coadd_xy = np.median(sigma_xys)

sigma_multi_flux = np.median(sigma_multi_fluxes)
sigma_multi_xy = np.median(sigma_multi_xys)

if args.template:
    # use errors from noisy measurements
    sigma_multi_flux *= args.factor
    sigma_multi_xy *= args.factor
    table_multi = bfd.TemplateTable(args.weight_n, args.weight_sigma, args.sn_min, sigma_multi_xy, sigma_multi_flux,
                                    args.sigma_step, args.sigma_max)
    sigma_coadd_flux *= args.factor
    sigma_coadd_xy *= args.factor
    table_coadd = bfd.TemplateTable(args.weight_n, args.weight_sigma, args.sn_min, sigma_coadd_xy, sigma_coadd_flux,
                                    args.sigma_step, args.sigma_max)

else:
    table_multi = bfd.TargetTable(args.weight_n, args.weight_sigma, cov=None)
    table_coadd = bfd.TargetTable(args.weight_n, args.weight_sigma, cov=None)


for i in range(args.ngal):
    if i%(args.ngal/10)==0 and i>0:
        print("%d%% done"% int(100.0*i/args.ngal))

    obs_list = sims(psf_seed=args.psf_seed)
    coadd_image = coaddsim.CoaddImages(obs_list, interp=args.interp, flat_wcs=args.flat_wcs, noise_n=args.noise_n)
    coadd = coadd_image.get_mean_coadd()

    bfd_multi = BfdObs(obs_list, weight, id=i, nda=1./args.ngal)
    bfd_coadd = BfdObs(coadd,    weight, id=i, nda=1./args.ngal, compute_noise_ps=args.use_noise_ps)

    if args.template:
        templates = bfd_multi.moment.make_templates(sigma_multi_xy, sn_min=args.sn_min,
                                                    sigma_flux=sigma_multi_flux, sigma_step=args.sigma_step,
                                                    sigma_max=args.sigma_max, xy_max=args.xy_max)
        for tmpl in templates:
            if tmpl is not None:
                table_multi.add(tmpl)

        templates = bfd_coadd.moment.make_templates(sigma_coadd_xy, sn_min=args.sn_min,
                                                    sigma_flux=sigma_coadd_flux, sigma_step=args.sigma_step,
                                                    sigma_max=args.sigma_max, xy_max=args.xy_max)
        for tmpl in templates:
            if tmpl is not None:
                table_coadd.add(tmpl)
    else:
        xyshift, error, msg = bfd_multi.moment.recenter()
        if error:
            #print 'Error:',i,msg,xyshift
            table_multi.addLost()
        else:
            table_multi.add(bfd_multi.moment.get_moment(0,0), xy=xyshift, id=i,
                            covgal=bfd_multi.moment.get_covariance())

        xyshift, error, msg = bfd_coadd.moment.recenter()
        if error:
            #print 'Error:',i,msg,xyshift
            table_coadd.addLost()
        else:
            table_coadd.add(bfd_coadd.moment.get_moment(0,0), xy=xyshift, id=i,
                            covgal=bfd_coadd.moment.get_covariance())

if args.template:
    table_multi.save('%s/%s_template_multi%s.fits'%(args.output_dir,args.name,args.label))
    table_coadd.save('%s/%s_template_coadd%s.fits'%(args.output_dir,args.name,args.label))
else:
    table_multi.save('%s/%s_multi%s.fits'%(args.output_dir,args.name,args.label))
    table_coadd.save('%s/%s_coadd%s.fits'%(args.output_dir,args.name,args.label))


