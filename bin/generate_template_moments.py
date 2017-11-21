#! /usr/bin/env python

import nsim
import galsim
import bfd
import numpy as np
import coaddsim
import argparse
from bfd_coadd import BfdObs
from astropy.table import Table


parser = argparse.ArgumentParser(description='Run single file')
parser.add_argument('--weight_n',default=4,type=int, help='nfor weight function')
parser.add_argument('--weight_sigma',default=1.,type=np.float64, help='sigma for weight function')
parser.add_argument('--sigma_step',default=1.2,type=float, help='sigma step for prior')
parser.add_argument('--sigma_max',default=5.5,type=float, help='sigma max')
parser.add_argument('--xy_max',default=4,type=float, help='xy max')
parser.add_argument('--sn_min',default=5,type=float, help='min s/n')
parser.add_argument('--ngal',default=1000,type=int, help='sigma max')
parser.add_argument('--noise_range', help='Noise is using range, needed for templates', action='store_const', const=True,default=None)
parser.add_argument('--file',default='bdk',type=str, help='sigma max')
parser.add_argument('--name',default='bdk',type=str, help='sigma max')
parser.add_argument('--label',default='',type=str, help='label')
parser.add_argument('--factor',default=100,type=float, help='noise factor for template')
parser.add_argument('--sigma',default=None,type=float, help='Use this sigma instead of value from file')
parser.add_argument('--output_dir',default='.', help='output directory')
parser.add_argument('--use_noise_ps', dest='use_noise_ps', default=False, action='store_true')
parser.add_argument('--seed',default=None,type=int, help='use this seed')
parser.add_argument('--psf_seed',default=-1,type=int, help='use this seed')
parser.add_argument('--coadd_cov_file',default=None,help='use covariance file')
parser.add_argument('--multi_cov_file',default=None,help='use covariance file')
parser.add_argument('--interp',default='lanczos3',type=str, help='interpolation kernel')
parser.add_argument('--flat_wcs',dest='flat_wcs', default=False, action='store_true')
args = parser.parse_args()


if args.seed is None:
    seed = np.random.randint(0,2**30)
else:
    seed = args.seed
weight = bfd.KSigmaWeight(args.weight_n, args.weight_sigma)

table_multi = None
table_coadd = None

sims = nsim.sime.Sim('%s'%(args.file), seed)

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
print new_config

new_config['shear'] = [0.0, 0.0]
sims = nsim.sime.Sim(new_config,seed)

coadd_cov = Table.read(args.coadd_cov_file)
multi_cov = Table.read(args.multi_cov_file)


tables_coadd=[]
sigmas_coadd_xy=[]
sigmas_coadd_flux=[]
coadd_bins=[]
for i in range(len(coadd_cov)):
    sigma_coadd_flux = np.sqrt(coadd_cov['cov_even'][i][0])
    sigma_coadd_xy = np.sqrt(coadd_cov['cov_odd'][i][0])
    table_coadd = bfd.TemplateTable(args.weight_n, args.weight_sigma, args.sn_min, sigma_coadd_xy, sigma_coadd_flux,
                                    args.sigma_step, args.sigma_max)
    tables_coadd.append(table_coadd)
    sigmas_coadd_xy.append(sigma_coadd_xy)
    sigmas_coadd_flux.append(sigma_coadd_flux)
    coadd_bins.append(coadd_cov['bin'][i])

tables_multi=[]
sigmas_multi_xy=[]
sigmas_multi_flux=[]
multi_bins=[]
for i in range(len(multi_cov)):
    sigma_multi_flux = np.sqrt(multi_cov['cov_even'][i][0])
    sigma_multi_xy = np.sqrt(multi_cov['cov_odd'][i][0])
    table_multi = bfd.TemplateTable(args.weight_n, args.]weight_sigma, args.sn_min, sigma_multi_xy, sigma_multi_flux,
                                    args.sigma_step, args.sigma_max)
    tables_multi.append(table_multi)
    sigmas_multi_xy.append(sigma_multi_xy)
    sigmas_multi_flux.append(sigma_multi_flux)
    multi_bins.append(multi_cov['bin'][i])



for i in range(args.ngal):
    if i%(args.ngal/10)==0 and i>0:
        print "%d%% done"% int(100.0*i/args.ngal)

    obs_list = sims(psf_seed=args.psf_seed)
    coadd_image = coaddsim.CoaddImages(obs_list, interp=args.interp, flat_wcs=args.flat_wcs)
    coadd = coadd_image.get_mean_coadd()

    bfd_multi = BfdObs(obs_list, weight, id=i, nda=1./args.ngal)
    bfd_coadd = BfdObs(coadd,    weight, id=i, nda=1./args.ngal, compute_noise_ps=args.use_noise_ps)

    for table,sigma_xy,sigma_flux in zip(tables_coadd, sigmas_coadd_xy, sigmas_coadd_flux):
        templates = bfd_coadd.moment.make_templates(sigma_xy, sn_min=args.sn_min,
                                                    sigma_flux=sigma_flux, sigma_step=args.sigma_step,
                                                    sigma_max=args.sigma_max, xy_max=args.xy_max)
        for tmpl in templates:
            if tmpl is not None:
                table.add(tmpl)

    for table,sigma_xy,sigma_flux in zip(tables_multi, sigmas_multi_xy, sigmas_multi_flux):
        templates = bfd_multi.moment.make_templates(sigma_xy, sn_min=args.sn_min,
                                                    sigma_flux=sigma_flux, sigma_step=args.sigma_step,
                                                    sigma_max=args.sigma_max, xy_max=args.xy_max)
        for tmpl in templates:
            if tmpl is not None:
                table.add(tmpl)


for bin,table in zip(coadd_bins,tables_coadd):
    table.save('%s/%s_template_coadd%s_%s.fits'%(args.output_dir,args.name,args.label,bin))

for bin,table in zip(multi_bins,tables_multi):
    table.save('%s/%s_template_multi%s_%s.fits'%(args.output_dir,args.name,args.label,bin))





