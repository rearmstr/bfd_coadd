import bfd
import numpy as np


class BfdObs(object):
    '''Construct a BFD MultiObject moment object from a list of nsim observations'''
    def __init__(self, obs_list, wt, id=0, nda=1.0):
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
            noise = np.sqrt(1./(np.mean(obs.weight)))
            kdata.append(bfd.simpleImage(obs.image, origin, obs.psf.image, wcs=wcs, pixel_noise=noise))

        self.moment = bfd.MultiMomentCalculator(kdata, wt, id=id, nda=nda)
