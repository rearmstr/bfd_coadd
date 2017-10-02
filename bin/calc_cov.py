import numpy as np
from astropy.table import Table
import glob
import argparse

parser = argparse.ArgumentParser(description='Run program')
parser.add_argument('--files',default="",help='files')
parser.add_argument('--var',default=0.1,type=float,help='allowed variation')
parser.add_argument('--low',default=1,type=float,help='minimnum percentile')
parser.add_argument('--high',default=99,type=float,help='minimnum percentile')
parser.add_argument('--outfile',default="out.fits",help='file to write out')
args = parser.parse_args()


files = glob.glob(args.files)
cov_even = []
cov_odd = []

files_read = 0
for file in files:
    try:
        f=Table.read(file)
    except Exception as e:
        print e
        continue

    cov_even.extend(f['cov_even'])
    cov_odd.extend(f['cov_odd'])
    files_read +=1

print 'read',files_read,'files'

cov_even = np.array(cov_even)
cov_odd = np.array(cov_odd)

lmin = np.percentile(np.sqrt(cov_even[:,0]),args.low)
lmax = np.percentile(np.sqrt(cov_even[:,0]),args.high)
med = np.median(np.sqrt(cov_even[:,0]))

diff = (lmax-lmin)/med
bins = int(np.ceil(diff/(args.var)))
binsize = (lmax-lmin)/bins
final_bins=[]
final_ave_even=np.zeros((bins,len(cov_even[0])))
final_ave_odd=np.zeros((bins,len(cov_odd[0])))

flux_cov_min = np.zeros(bins)
flux_cov_max = np.zeros(bins)

print 'full range',lmin,lmax
for i in range(bins):
    cmin = i*binsize + lmin
    cmax = (i+1)*binsize + lmin
    print 'bin',i+1,'min',cmin,'max',cmax

    mask = (np.sqrt(cov_even[:,0])>cmin) & (np.sqrt(cov_even[:,0])<cmax)
    ave_even = np.mean(cov_even[mask],axis=0)
    ave_odd = np.mean(cov_odd[mask],axis=0)

    print ' even',ave_even
    print ' odd',ave_odd
    final_ave_even[i] = ave_even
    final_ave_odd[i] = ave_odd
    flux_cov_min[i] = cmin
    flux_cov_max[i] = cmax
    final_bins.append('bin%d'%(i+1))

cat = Table([final_ave_even,final_ave_odd,final_bins,flux_cov_min,flux_cov_max],
            names=['cov_even','cov_odd','bin','flux_cov_min','flux_cov_max'])

cat.write(args.outfile,format='fits',overwrite=True)
