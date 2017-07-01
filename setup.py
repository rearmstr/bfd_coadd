from distutils.core import setup


setup(
    name="bfd_coadd",
    packages=['bfd_coadd'],
    version="0.1",
    scripts=['bfd_coadd/submit_moments.py', 'bfd_coadd/submit_lsst-dev_moments.py',
             'bfd_coadd/submit_moments.py', 'bfd_coadd/generate_moments.py',
             'bfd_coadd/generate_mp_moments.py'
            ]
)




