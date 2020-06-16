import astroquery
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from astroquery.vo_conesearch import ConeSearch
from matplotlib.patches import Circle, Rectangle


__all__ = ['GuideStars']


class GuideStars():
    """ Finding Guide Stars based on selected catalogues by astroquery.
    """
    def __init__(self, ra, dec, radius=0.1, unit='deg'):
        '''
        Parameters
        ----------
        ra, dec : float, str, `~astropy.Quantity`
            The RA and DEC position of the center of the FOV for
            cone-search. Currently only RA and DEC is supported for
            brevity of the API. If float, it must be in the unit of
            ``unit``.
        radius : float, `~astropy.Quantity`, optional
            The cone search radius. If float, it must be in the unit of
            ``unit``.
        unit : str, `~astropy.unit`, optional.
            The unit for ``ra``, ``dec``, and ``radius`` if they're
            floats.

        Attributes
        ----------
        center : The center of FOV in `~astropy.coordinates.SkyCoord`.
        radius : The search radius
        center_frame : the `~astropy.coordinates.SkyOffsetFrame`
        centered at the ``self.center``.
        '''
        self.center = SkyCoord(ra=ra, dec=dec, unit=unit)
        if isinstance(radius, u.Quantity):
            self.radius = radius.to(unit)
        else:
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
        ''' Queries and organizes the query result
        Parameters
        ----------
        mag_cut : float, optional.
            The maximum magnitude (in ``'Mag'`` column of the queried
            result in pure astroquery) to use as guide star. If
            ``None``, no such cut will be applied.

        Note
        ----
        Internally it will add the following columns:
            * ``'_r'``: Angular distance from the center
            * ``'offset_x'``, ``'offset_y'``: angular offsets with
            respect to ``self.center`` in the ``self.center_frame``
            frame. This corrects the high-DEC effect (at DEC = 89 deg, 1
            deg separation changes RA value dramatically by ~ 1/cos(DEC)
            factor.)
        ``self.mag_cut`` is the ``mag_cut`` given, and
        ``self.mag_ref`` is that ``mag_cut``. If ``mag_cut`` is
        ``None``, ``self.mag_ref`` is the maximum value of the
        magnitudes of the guide stars.
        '''
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
        self.mag_cut = mag_cut
        self.mag_ref = mag_cut if mag_cut is not None else np.max(mags)

    @staticmethod
    def m2s(mag, mag_ref, size_factor=1.):
        ''' Change magnitude to matplotlib scatter marker size.
        Parameters
        ----------
        mag_ref : float
            The reference magnitude which marker size is
            ``size_factor``.
        size_factor : float, optional.
            The factor which multicatively scales the marker size.
        '''
        return size_factor*10**(-0.4*(mag - mag_ref))

    @staticmethod
    def s2m(size, mag_ref, size_factor=1.):
        ''' Change matplotlib scatter marker size to magnitude.
        Parameters
        ----------
        mag_ref : float
            The reference magnitude which marker size is
            ``size_factor``.
        size_factor : float, optional.
            The factor which multicatively scales the marker size.
        '''
        return mag_ref - 2.5*np.log10(size/size_factor)

    def quickplot(
        self, ax, marker='o', color='k', size_factor=20., size=None,
        alpha=0.3, scat_kw={},
        cross=True, cross_kw={'color': 'k', 'lw': 0.5},
        unit='arcmin', invert_lon=True, query_circle=True,
        circ_kw={'facecolor': 'none', 'edgecolor': 'k', 'lw': 1, 'ls': ':'},
        fov_size=None, fov_angle=0,
        fov_kw={'facecolor': 'none', 'edgecolor': 'k', 'lw': 1},
        rect_kw={'facecolor': 'none', 'edgecolor': 'k', 'lw': 1},
        num_show_mag=5, mag_kw={'fontsize': 10},
        legend=True, legend_kw={'loc': 1, 'fontsize': 9}
    ):
        ''' Quickly draw FoV plot. Maybe no need for documentation.
        '''
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

        rad = self.radius.to_value(unit)

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

        if query_circle:
            circ = Circle((0, 0), rad, **circ_kw)
            ax.add_patch(circ)

        width = 0.
        height = 0.
        if fov_size is not None:
            fov_size = np.atleast_1d(fov_size)
            if fov_size.ndim != 1:
                raise ValueError("fov_size must be scalar or 1-D")
            if fov_size.shape[0] == 1:
                width, height = fov_size[0], fov_size[0]
            elif fov_size.shape[0] == 2:
                width, height = fov_size[0], fov_size[1]
            else:
                raise ValueError("fov_size must be size of 1 or 2.")

            rect = Rectangle((-width/2, -height/2), width=width, height=height,
                             angle=fov_angle, **fov_kw)
            ax.add_patch(rect)

        # iterate only first ``num_show_mag`` rows
        if show_mag:
            for i in range(num_show_mag):
                ax.text(offsets_x[i], offsets_y[i],
                        s=f"{mags[i]:.2f}", **mag_kw)

        ax.set(
            xlim=np.array([-1., 1.])*max(width/2., rad),
            ylim=np.array([-1., 1.])*max(height/2., rad),
        )

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

        ax.set_aspect('equal')
        print(f"FoV = ({width:.1f}, {height:.1f}) {unit}")
