from distutils.core import setup


setup(
    name="bfd_coadd",
    packages=['bfd_coadd'],
    version="0.1",
    scripts=['bin/submit_moments.py', 'bin/submit_lsst-dev_moments.py',
             'bin/submit_integrate.py', 'bin/generate_moments.py',
             'bin/generate_mp_moments.py', 'bin/submit_integrateN.py'
            ]
)




