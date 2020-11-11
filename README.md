# ysobsplanpy

A simple small package to make a finding chart of your object using guide star catalogs.

## Requirements
```
numpy
astropy
astroquery (better to use the latest version)
matplotlib (optional but recommended)
```

## Usage

```python
%config InlineBackend.figure_format = 'retina'
from matplotlib import pyplot as plt
from matplotlib import rcParams
from ysobsplanpy import GuideStars

# We need to do it in a separate cell. See:
# https://github.com/jupyter/notebook/issues/3385
plt.style.use('default')
rcParams.update({'font.size':12})

# a random object SDSS J160336.42+155906.8 -- Quasar
g = GuideStars(ra=240.90177, dec=15.985262, radius=0.3, unit='deg')
g.query(mag_cut=13)
fig, axs = plt.subplots(1, 1, figsize=(7, 7),
                        sharex=False, sharey=False, gridspec_kw=None)

# Try g.quickplot? and change parameters as you like.
g.quickplot(axs, num_show_mag=5, fov_size=21, unit='arcmin')

plt.tight_layout()
fig.align_ylabels(axs)
fig.align_xlabels(axs)
plt.show()
```



![](figs/example.png)


or

```python
%config InlineBackend.figure_format = 'retina'
from matplotlib import pyplot as plt
from matplotlib import rcParams
import ysobsplanpy as op
# We need to do it in a separate cell. See:
# https://github.com/jupyter/notebook/issues/3385
plt.style.use('default')
rcParams.update({'font.size':12})
pos = dict(
    ceres1=dict(ra="22:25:01.46", dec="-25:11:41.2"),
    ceres2=dict(ra="22:25:11.30", dec="-24:34:44.3"),
#     ceres2=dict(ra="22:25:11.30", dec="-24:54:44.3"),
    pallas1=dict(ra="18:56:42.74", dec="+05:06:22.8"),
#     pallas2=dict(ra="20:06:42.28", dec="-00:27:19.8"),
    suhyun=dict(ra="01:17:33", dec="-04:44:44"),
)

fig, axs = plt.subplots(2, 2, figsize=(12, 12),
                        sharex=False, sharey=False, gridspec_kw=None)

for i, (name, v) in enumerate(pos.items()):
    ax = axs.flatten()[i]
    gs = op.GuideStars(ra=v['ra'], dec=v['dec'], radius=0.5)
    gs.query(mag_cut=20)
    gs.quickplot(ax)
    ax.set_title(name)


plt.tight_layout()
fig.align_ylabels(axs)
fig.align_xlabels(axs)
plt.show()

# plt.savefig(figoutpath, dpi=300, bbox_inches = "tight")
```
