import astroquery
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from astroquery.vo_conesearch import ConeSearch
from matplotlib.patches import Circle


__all__ = ['GuideStars']


class GuideStars():
    def __init__(self, ra, dec, radius=0.1, unit='deg'):
        self.center = SkyCoord(ra=ra, dec=dec, unit=unit)
        self.radius = radius*u.Unit(unit)
        self.center_frame = SkyOffsetFrame(origin=self.center)

    def __str__(self):
        expl = (
            "Cone-search for guide stars using astroquery.vo_conesearch "
            + f"(astroquery v. {astroquery.__version__})."
            + f" Search radius {self.radius} around {self.center}"
        )
        return expl

    def query(self, mag_cut=13.):
        self.queried = ConeSearch.query_region(self.center, radius=self.radius)
        if mag_cut is not None:
            magmask = self.queried['Mag'] > mag_cut
            self.queried = self.queried[~magmask]
        # self.queried['dra'] = self.queried['ra'].to(u.deg) - self.center
        self.stars = SkyCoord(self.queried['ra'], self.queried['dec'])

        self.queried['_r'] = self.center.separation(self.stars)

        dist_offset = self.stars.transform_to(self.center_frame)
        self.queried['offset_x'] = dist_offset.lon
        self.queried['offset_y'] = dist_offset.lat
        mags = self.queried['Mag']
        self.magcut = mag_cut
        self.mag_ref = mag_cut if mag_cut is not None else np.min(mags)

    @staticmethod
    def m2s(mag, mag_ref, size_factor=1.):
        return size_factor*10**(-0.4*(mag - mag_ref))

    @staticmethod
    def s2m(size, mag_ref, size_factor=1.):
        return mag_ref - 2.5*np.log10(size/size_factor)

    def quickplot(
            self, ax, marker='o', color='k', size_factor=20., size=None,
            alpha=0.3, scat_kw={},
            cross=True, cross_kw={'color': 'k', 'lw': 0.5},
            unit='arcmin', invert_lon=True, fov_circle=True,
            fov_kw={'facecolor': 'none', 'edgecolor': 'k', 'lw': 1},
            num_show_mag=5, mag_kw={'fontsize': 10},
            legend=True, legend_kw={'loc': 1, 'fontsize': 9}
    ):
        show_mag = False
        if num_show_mag is not None:
            show_mag = True
            n_queried = len(self.queried)
            if num_show_mag < 0 or num_show_mag > n_queried:
                num_show_mag = n_queried
            self.queried.sort('Mag')

        if size is None:
            size = self.m2s(self.queried['Mag'].copy(), self.mag_ref,
                            size_factor=size_factor)
        else:
            size_factor = 1.

        offsets_x = np.array(self.queried['offset_x'].to(unit))
        offsets_y = np.array(self.queried['offset_y'].to(unit))
        mags = self.queried['Mag']
        ax.scatter(
            offsets_x,
            offsets_y,
            marker=marker,
            color=color,
            s=size,
            alpha=alpha,
            **scat_kw
        )
        ax.scatter(0, 0, marker='s', color='r', facecolors='none',
                   s=size_factor)
        ax.set(
            xlabel=f"Projected offset (longitude; {unit})",
            ylabel=f"Projected offset (latitude; {unit})"
        )

        if cross:
            ax.axhline(0, **cross_kw)
            ax.axvline(0, **cross_kw)

        if fov_circle:
            rad = self.radius.to_value(unit)
            circ = Circle((0, 0), rad, **fov_kw)
            ax.add_patch(circ)
            ax.set(
                xlim=(-rad, +rad),
                ylim=(-rad, +rad)
            )

        # iterate only first ``num_show_mag`` rows
        if show_mag:
            for i in range(num_show_mag):
                ax.text(offsets_x[i], offsets_y[i],
                        s=f"{mags[i]:.2f}", **mag_kw)

        if invert_lon:
            ax.set(
                xlim=(ax.get_xlim()[::-1])
            )

        if legend:
            mm = np.arange(np.floor(np.min(mags)), self.mag_ref + 0.1, 1.)
            # plt.plot has different size notation... F...k....
            ss = np.sqrt(self.m2s(mm, self.mag_ref, size_factor))
            for m, s in zip(mm, ss):
                ax.plot(np.nan, np.nan, marker=marker, ls='', alpha=alpha,
                        color=color, label=m, ms=s)
            ax.legend(title="Mag", **legend_kw)
