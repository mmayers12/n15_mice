#! /usr/bin/env python3

"""
Generate a xcorr vs deltcn plot using all sqt files in the current dir

"""

import glob, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from itertools import chain
from . import blazmass_tools
import numpy as np
import sys

if len(sys.argv) == 2:
    sqt_path = os.path.expanduser(sys.argv[1])
else:
    sqt_path = os.getcwd()

sqt_files = glob.glob(os.path.join(sqt_path, '*.sqt'))
names = set(filename.rsplit('.')[0].rsplit('_',1)[0] for filename in sqt_files)
if len(names) == 1:
    name = os.path.basename(list(names)[0])
else:
    name = os.path.basename(sqt_path)

sqt_chunker = chain(*[blazmass_tools.sqt_chunker(sqt_file) for sqt_file in sqt_files])
xcorr, deltcn, reverse = zip(*[(x['Xcorr'],x['deltCN'],x['matches'][0]['reverse']) for x in sqt_chunker])
reverse = np.array(reverse)
deltcn = np.array(deltcn)
xcorr = np.array(xcorr)

#%%
colors = {}
colors['forward_human'] = '#468966'
colors['forward_nonhuman'] = '#FFF0A5'
colors['reverse'] = '#FFB03B'
colors['filtered'] = '#B64926'
circlearea = 4
legend_marker_size = 5
opacity = 0.8
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
ax.scatter(xcorr[~reverse], deltcn[~reverse], c = 'r', alpha = opacity, label = 'Forward PSMs', lw=0.1, s=circlearea)
ax.scatter(xcorr[reverse], deltcn[reverse], c = 'b', alpha = opacity, label = 'Reverse PSMs', lw=0.1, s=circlearea)
ax.set_xlabel('XCorr',fontsize=24)
ax.set_ylabel('DeltCN',fontsize=24)
ax.set_xlim(left = 0, right = 10)
ax.set_ylim(bottom = 0, top = 1.0)
ax.legend(loc='upper right', markerscale=4, fancybox=True)
plt.title("XCorr vs. DeltCN\n{}".format(name))
plt.tight_layout()
plt.savefig(os.path.join(sqt_path,"xcorr_deltcn_{}.png".format(name)), orientation='landscape', dpi=300)
plt.savefig(os.path.join(sqt_path,"xcorr_deltcn_{}.pdf".format(name)), orientation='landscape')
