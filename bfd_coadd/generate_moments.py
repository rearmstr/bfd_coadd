import nsim
import galsim
import bfd
import numpy as np
import coaddsim
import argparse


class BfdObs(object):
    '''Construct a BFD MultiObject moment object from a list of nsim observations'''
    def __init__(self, obs_list, wt, id=0, nda=1.0, use_offset=True):
        if isinstance(obs_list, list) is False:
            obs_list = [obs_list]
        kdata = []
        for obs in obs_list:

            jacobian = np.array([ [obs.jacobian.dudcol, obs.jacobian.dudrow],
                                  [obs.jacobian.dvdcol,obs.jacobian.dvdrow] ])
            center = (np.array(obs.image.shape)-1.0)/2.0
            yxref =  center
            if use_offset:
                center += obs.meta['offset_pixels']

            xyref = [yxref[1], yxref[0]]
            uvref = (0., 0.)
            wcs = bfd.WCS(jacobian, xyref=xyref, uvref=uvref)
            origin = uvref
            noise = np.sqrt(1./(np.mean(obs.weight)))
            kdata.append(bfd.simpleImage(obs.image, origin, obs.psf.image, wcs=wcs, pixel_noise=noise))

        self.moment = bfd.MultiMomentCalculator(kdata, wt, id=id, nda=nda)

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
parser.add_argument('--label',default='',type=str, help='label')
parser.add_argument('--factor',default=100,type=float, help='noise factor for template')
parser.add_argument('--output_dir',default='.', help='output directory')

args = parser.parse_args()



seed=int(np.random.rand()*100000000)
weight = bfd.KSigmaWeight(args.ngal,args.weight_sigma)

table_multi = None
table_coadd = None

sims = nsim.sime.Sim('%s'%(args.file), seed)

# Set noise higher and shear to zero and recreate sims
if args.template:
    new_config={}
    for key,value in sims.items():
        new_config[key] = value

    new_config['images']['noise'] = sims['images']['noise']/args.factor
    new_config['shear'] = [0.0, 0.0]
    sims = nsim.sime.Sim(new_config,seed)

obs_test = sims()
bfd_test = BfdObs(obs_test, weight, id=0, nda=1./args.ngal)
cov_test = bfd_test.moment.get_covariance()

sigma_flux = np.sqrt(cov_test[0][0,0])
sigma_xy = np.sqrt(cov_test[1][0,0])

if args.template:
    # use errors from noisy measurements
    sigma_flux *= args.factor
    sigma_xy *= args.factor
    table_multi = bfd.TemplateTable(args.weight_n, args.weight_sigma, args.sn_min, sigma_xy, sigma_flux,
                                    args.sigma_step, args.sigma_max)
    table_coadd = bfd.TemplateTable(args.weight_n, args.weight_sigma, args.sn_min, sigma_xy, sigma_flux,
                                    args.sigma_step, args.sigma_max)

else:
    table_multi = bfd.TargetTable(args.weight_n, args.weight_sigma, cov=None)
    table_coadd = bfd.TargetTable(args.weight_n, args.weight_sigma, cov=None)


for i in range(args.ngal):
    if i%(args.ngal/10)==0 and i>0:
        print "%d%% done"% int(100.0*i/args.ngal)

    obs_list = sims()
    coadd_image = coaddsim.CoaddImages(obs_list)
    coadd = coadd_image.get_mean_coadd(False)

    bfd_multi = BfdObs(obs_list, weight, id=i, nda=1./args.ngal)
    bfd_coadd = BfdObs(coadd,    weight, id=i, nda=1./args.ngal, use_offset=False)

    if args.template:
        templates = bfd_multi.moment.make_templates(sigma_xy, sn_min=args.sn_min,
                                                    sigma_flux=sigma_xy, sigma_step=args.sigma_step,
                                                    sigma_max=args.sigma_max, xy_max=args.xy_max)
        for tmpl in templates:
            if tmpl is not None:
                table_multi.add(tmpl)

        templates = bfd_coadd.moment.make_templates(sigma_xy, sn_min=args.sn_min,
                                                    sigma_flux=sigma_flux, sigma_step=args.sigma_step,
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


#nohup python submit_integrate.py --njobs 50 --nfiles 1000 --dir target_1p/ --template target_1p_t10000_template_multi_0.fits --template_label t10000_multi  --label target_1p_multi --config test.param  --hours 8 --max_jobs 15 --mins 3 >&log.target_1p.multi&

#nohup python submit_integrate.py --njobs 50 --nfiles 1000 --dir target_1p/ --template target_1p_t10000_template_coadd_0.fits --template_label t10000_coadd  --label target_1p_coadd --config test.param  --hours 8 --max_jobs 15 --mins 3 >&log.target_1p.coadd&
