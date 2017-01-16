 # -*- coding: utf-8 -*-
from __future__ import print_function
import matplotlib.pylab as plt
import numpy as np
import matplotlib as mpl
import matplotlib.lines as mlines
from pylab import Rectangle
import os
from obspy import read
mpl.use('Agg')
dists = np.load('ds.npy')
mags = np.load('mags.npy')
dist_deg = []

for d in range(0,len(dists)):
	dist_deg.append(dists[d]/111.11)



fig = plt.figure(figsize=(20/2.54,7/2.54))
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10

### VORSCHLAG 1

bins=[0,0,0,0,0]

for m in range(0,len(mags)):
	if ((mags[m] >= 4) and (mags[m] < 5)):
		bins[0] += 1
	elif ((mags[m] >= 5) and (mags[m] < 6)):
		bins[1] += 1
	elif ((mags[m] >= 6) and (mags[m] < 7)):
		bins[2] += 1
	elif ((mags[m] >= 7) and (mags[m] < 8)):
		bins[3] += 1
	elif (mags[m] >=8):
		bins[4] += 1
print(bins)

ax1 = plt.subplot2grid((1, 10), (0, 0), colspan=7)
ax1.axhline(4, linewidth=.3, color='k')
ax1.axhline(5, linewidth=.3, color='k')
ax1.axhline(6, linewidth=.3, color='k')
ax1.axhline(7, linewidth=.3, color='k')
ax1.axhline(8, linewidth=.3, color='k')
ax1.scatter(dist_deg, mags, s=.1)
ax1.set_xlabel(u'epicentral distance [°]', fontsize=11, labelpad=3)
ax1.set_ylabel('moment magnitude', fontsize=11)
ax1.set_xlim(-40,180)
ax1.set_ylim(3,9.2)
ax1.set_yticks([4,5,6,7,8,9])
ax1.set_xticks([0,30,60,90,120,150])
ax1.text(s='# events', x=-35, y=9.0, fontsize=9)
ax1.text(s=str(bins[0]), x=-35, y=4.35, fontsize=9)
ax1.text(s=str(bins[1]), x=-35, y=5.35, fontsize=9)
ax1.text(s=str(bins[2]), x=-35, y=6.35, fontsize=9)
ax1.text(s=str(bins[3]), x=-35, y=7.35, fontsize=9)
ax1.text(s=str(bins[4]), x=-35, y=8.35, fontsize=9)


### Vorschlag 2

bins2 = np.zeros(53, dtype=np.int)
min=4.0
max=4.1

for counter in range(0,len(bins2)):
	for m2 in range(0,len(mags)):
		if ((mags[m2] >= min) and (mags[m2] < max)):
			bins2[counter] += 1
	min+=0.1
	max+=0.1

ax2 = plt.subplot2grid((1, 10), (0, 7), colspan=3, sharey=ax1)
mag_bins = np.linspace(4.0,9.2,53)
ax2.barh(bottom=mag_bins-.05, width=bins2, height=0.11, color='k', edgecolor='w')
ax2.text(s='total: '+ str(len(mags)), x=200, y=9.0, fontsize=9)
ax2.set_xscale('log')
ax2.set_xlim(0.6,11000)
ax2.set_ylim(3.5,9.5)
ax2.yaxis.tick_right()
#ax2.set_yticks([])
#ax2.set_xticks([4,5,6,7,8,9])
#ax2.text(s='events:\nJUL 2007 - JUL 2016\ntotal number:\n'+str(len(mags)), x=6.5, y=500, fontsize=9)
ax2.set_xlabel('# events', fontsize=11, labelpad=1)
#ax2.set_ylabel('# processed events', fontsize=11)
ax2.axhline(4, linewidth=.3, color='k')
ax2.axhline(5, linewidth=.3, color='k')
ax2.axhline(6, linewidth=.3, color='k')
ax2.axhline(7, linewidth=.3, color='k')
ax2.axhline(8, linewidth=.3, color='k')
plt.subplots_adjust(right=0.96, wspace=0.12,bottom=0.18, left=0.07, top=0.9)
#plt.tight_layout()
plt.show()
#plt.savefig('Latex_pp/dist_mag_merged.pdf')