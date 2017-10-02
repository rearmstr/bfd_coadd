import bfd
import numpy as np

def calc_noise_ps(image, min_k=0, max_k=4.4, dk=0.6):
    k_xaxis = np.fft.fftfreq(image.shape[1])*2.*np.pi
    k_yaxis = np.fft.fftfreq(image.shape[0])*2.*np.pi
    kx, ky = np.meshgrid(k_xaxis, k_yaxis)

    absk = np.sqrt(kx**2 + ky**2)
    ps = np.abs(np.fft.fft2(image))**2

    kbins = np.arange(min_k, max_k, dk)
    ave_k=np.zeros(len(kbins)-1)
    ave_ps=np.zeros(len(kbins)-1)
    ave_ps_err=np.zeros(len(kbins)-1)

    for i in range(len(kbins)-1):
        min_k, max_k = kbins[i], kbins[i+1]
        mask = (absk >= min_k) & (absk < max_k)
        ave_ps[i] = np.mean(ps[mask])
        ave_ps_err[i] = np.std(ps[mask])/np.sqrt(np.sum(mask))
        ave_k[i] = (min_k+max_k)/2.
    ave_ps = ave_ps/(image.shape[0]*image.shape[1])

    return np.array(ave_k), ave_ps

class Interp:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __call__(self, v):
        return np.interp(v, self.x, self.y)

class BfdObs(object):
    '''Construct a BFD MultiObject moment object from a list of nsim observations'''
    def __init__(self, obs_list, wt, id=0, nda=1.0, compute_noise_ps=False):
        if isinstance(obs_list, list) is False:
            obs_list = [obs_list]
        kdata = []
        for obs in obs_list:

            jacobian = np.array([ [obs.jacobian.dudcol, obs.jacobian.dudrow],
                                  [obs.jacobian.dvdcol,obs.jacobian.dvdrow] ])
            xyref = (obs.jacobian.col0, obs.jacobian.row0)
            uvref = (0., 0.)
            wcs = bfd.WCS(jacobian, xyref=xyref, uvref=uvref)
            origin = uvref
            if compute_noise_ps is False:
                noise = np.sqrt(1./(np.mean(obs.weight)))
                kdata.append(bfd.simpleImage(obs.image, origin, obs.psf.image, wcs=wcs, pixel_noise=noise))
            else:
                k, ps = calc_noise_ps(obs.noise)
                interp = Interp(k, ps)
                kdata.append(bfd.simpleImage(obs.image, origin, obs.psf.image, wcs=wcs, noise_ps=interp))

        self.moment = bfd.MultiMomentCalculator(kdata, wt, id=id, nda=nda)
