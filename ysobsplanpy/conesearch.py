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

    def query(self, magcut=13.):
        self.queried = ConeSearch.query_region(self.center, radius=self.radius)
        if magcut is not None:
            magmask = self.queried['Mag'] > magcut
            self.queried = self.queried[~magmask]
        # self.queried['dra'] = self.queried['ra'].to(u.deg) - self.center
        self.stars = SkyCoord(self.queried['ra'], self.queried['dec'])
        dist_offset = self.stars.transform_to(self.center_frame)
        self.queried['offset_x'] = dist_offset.lon
        self.queried['offset_y'] = dist_offset.lat
        mags = self.queried['Mag']
        refmag = magcut if magcut is not None else np.min(mags)
        self.queried['brightness'] = 10**(-0.4*(mags - refmag))

    def quickplot(
            self, ax, marker='o', color='k', size_factor=20, size=None,
            alpha=0.3, cross=True, cross_kw={'color': 'k', 'lw': 0.5},
            unit='arcmin', invert_lon=True, fov_circle=True,
            fov_kw={'facecolor': 'none', 'edgecolor': 'k', 'lw': 1},
            num_show_mag=5, mag_kw={},
            **kwargs
    ):
        if size is None:
            size = size_factor*self.queried['brightness']

        if num_show_mag is not None:
            show_mag = True
            if num_show_mag < 0:
                num_show_mag = len(self.queried)
            self.queried.sort('Mag')

        offsets_x = np.array(self.queried['offset_x'].to(unit))
        offsets_y = np.array(self.queried['offset_y'].to(unit))
        ax.scatter(
            offsets_x,
            offsets_y,
            marker=marker,
            color=color,
            s=size,
            alpha=alpha,
            **kwargs
        )
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
                        s=f"{self.queried['Mag'][i]:.2f}", **mag_kw)

        if invert_lon:
            ax.set(
                xlim=(ax.get_xlim()[::-1])
            )
