import nsim
import galsim
import bfd
import numpy as np
import coaddsim
import argparse
import multiprocessing
import time


class BfdObs(object):
    '''Construct a BFD MultiObject moment object from a list of nsim observations'''
    def __init__(self, obs_list, wt, id=0, nda=1.0):
        if isinstance(obs_list, list) is False:
            obs_list = [obs_list]
        kdata = []
        for obs in obs_list:

            jacobian = np.array([ [obs.jacobian.dudcol, obs.jacobian.dudrow],
                                  [obs.jacobian.dvdcol,obs.jacobian.dvdrow] ])

            origin= (0., 0.)
            xyref = (obs.jacobian.row0, obs.jacobian.col0)
            uvref = (0., 0.)
            wcs = bfd.WCS(jacobian, xyref=xyref, uvref=uvref)
            noise = np.sqrt(1./(np.mean(obs.weight)))
            kdata.append(bfd.simpleImage(obs.image, origin, obs.psf.image, wcs=wcs, pixel_noise=noise))

        self.moment = bfd.MultiMomentCalculator(kdata, wt, id=id, nda=nda)



def worker(weight_n,weight_sigma,sigma_step,sigma_max,xy_max,sn_min,ngal,target,template,file,name,label,factor,output_dir,index):

    wname = multiprocessing.current_process().name
    seed=int(np.random.rand()*100000000)+index
    weight = bfd.KSigmaWeight(ngal,weight_sigma)

    table_multi = None
    table_coadd = None

    sims = nsim.sime.Sim('%s'%(file), seed)

    # Set noise higher and shear to zero and recreate sims
    if template:
        new_config={}
        for key,value in sims.items():
            new_config[key] = value

        new_config['images']['noise'] = sims['images']['noise']/factor
        new_config['shear'] = [0.0, 0.0]
        sims = nsim.sime.Sim(new_config,seed)

    obs_test = sims()
    bfd_test = BfdObs(obs_test, weight, id=0, nda=1./ngal)
    cov_test = bfd_test.moment.get_covariance()

    sigma_flux = np.sqrt(cov_test[0][0,0])
    sigma_xy = np.sqrt(cov_test[1][0,0])

    if template:
        # use errors from noisy measurements
        sigma_flux *= factor
        sigma_xy *= factor
        table_multi = bfd.TemplateTable(weight_n, weight_sigma, sn_min, sigma_xy, sigma_flux,
                                        sigma_step, sigma_max)
        table_coadd = bfd.TemplateTable(weight_n, weight_sigma, sn_min, sigma_xy, sigma_flux,
                                        sigma_step, sigma_max)

    else:
        table_multi = bfd.TargetTable(weight_n, weight_sigma, cov=None)
        table_coadd = bfd.TargetTable(weight_n, weight_sigma, cov=None)


    for i in range(ngal):
        if i%(ngal/10)==0 and i>0:
            print wname,"%d%% done"% int(100.0*i/ngal)
        obs_list = sims()
        coadd_image = coaddsim.CoaddImages(obs_list)
        coadd = coadd_image.get_mean_coadd(False)

        bfd_multi = BfdObs(obs_list, weight, id=i, nda=1./ngal)
        bfd_coadd = BfdObs(coadd, weight, id=i, nda=1./ngal)

        if template:
            templates = bfd_multi.moment.make_templates(sigma_xy, sn_min=sn_min,
                                                        sigma_flux=sigma_xy, sigma_step=sigma_step,
                                                        sigma_max=sigma_max, xy_max=xy_max)
            for tmpl in templates:
                if tmpl is not None:
                    table_multi.add(tmpl)

            templates = bfd_coadd.moment.make_templates(sigma_xy, sn_min=sn_min,
                                                        sigma_flux=sigma_flux, sigma_step=sigma_step,
                                                        sigma_max=sigma_max, xy_max=xy_max)
            for tmpl in templates:
                if tmpl is not None:
                    table_coadd.add(tmpl)
        else:
            xyshift, error, msg = bfd_multi.moment.recenter()
            if error:
                print 'Error:',i,msg,xyshift
                table_multi.addLost()
            else:
                table_multi.add(bfd_multi.moment.get_moment(0,0), xy=xyshift, id=i,
                                covgal=bfd_multi.moment.get_covariance())

            xyshift, error, msg = bfd_coadd.moment.recenter()
            if error:
                print 'Error:',i,msg,xyshift
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
    parser.add_argument('--output_dir',default='.', help='output directory')
    parser.add_argument('--start',default=0,type=int, help='output directory')
    parser.add_argument('--njobs',default=48,type=int, help='output directory')
    args = parser.parse_args()
    print args.start

    jobs = []
    for i in range(args.njobs):
        print 'starting',i
        label = '_%d'%(i+args.start)
        arg=(args.weight_n,args.weight_sigma,args.sigma_step,args.sigma_max,args.xy_max,args.sn_min,args.ngal,
              args.target,args.template,args.file,args.name,label,args.factor,args.output_dir,i+args.start)
        p = multiprocessing.Process(target=worker, args=arg)
        jobs.append(p)

    for j in jobs:
        j.start()

    for j in jobs:
        j.join()
