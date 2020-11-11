from astroquery.jplhorizons import Horizons, HorizonsClass
from astropy.time import Time
import numpy as np
from astroquery.jplsbdb import SBDB

__all__ = ["SmallBodyQuery"]


class SmallBodyQuery(HorizonsClass):
    def __init__(
        self, id, location=None, epochs=None, id_type='smallbody',
        start=None, stop=None, step=None, lon=None, lat=None, elevation=None
    ):
        if epochs is None:
            if (start is None) or (stop is None) or (step is None):
                raise ValueError(
                    "Either epochs or all three start, stop, step "
                    + "must be given.")
            epochs = dict(start=str(start), stop=str(stop), step=str(step))
        elif not isinstance(epochs, dict):
            raise ValueError(
                "This class is purely for observational strategy planning, so "
                + "the list-like epoch is not supported at the moment. "
                + "Use dict-like or start, stop, step."
            )

        if location is None:
            if (lon is None) or (lat is None) or (elevation is None):
                raise ValueError(
                    "Either location or all three lon[deg; west=negative], "
                    + "lat[deg; south=negative], levation[m] must be given.")
            location = dict(lon=float(lon), lat=float(lat),
                            elevation=float(elevation))

        super().__init__(id=id, location=location, epochs=epochs,
                         id_type=id_type)
        self.shortname = SBDB.query(self.id)['object']['shortname']

    def query(self, pandas=False, elmax_window=1, elmax_min_count=None,
              nightmask_el=10, *args, **kwargs):
        try:
            import bottleneck as bn
        except ImportError:
            raise ImportError("Bottleneck is required for this."
                              + "https://github.com/pydata/bottleneck")

        self.queried = self.ephemerides(*args, **kwargs)
        self.queried["ELmax"] = bn.move_max(self.queried['EL'],
                                            window=elmax_window,
                                            min_count=elmax_min_count)
        self.datetime = Time(self.queried['datetime_jd'], format='jd')
        self.plotdate = self.datetime.plot_date

        # angular speed = sqrt((v_RA*cosDEC)^2 + v_DEC^2) [arcsec/h]
        self.queried["ang_speed"] = np.sqrt(self.queried['RA_rate']**2
                                            + self.queried['DEC_rate']**2)
        if not pandas:
            self.queried["ang_speed"].unit = "arcsec/h"

        self.nightmask = ((self.queried['solar_presence'] == '')
                          & (self.queried['EL'] < nightmask_el))
        el_night = self.queried['EL']
        el_night[~self.nightmask] = np.nan
        self.queried["ELmax_night"] = bn.move_max(
            el_night,
            window=elmax_window,
            min_count=elmax_min_count
        )

    def plot_timeax(self, ax, toplot=["alpha"], plot_kw=None, **kwargs):
        try:
            import ysvisutilpy as yvu
        except ImportError:
            raise ImportError("This method requries ysvisutilpy at the moment."
                              + "https://github.com/ysBach/ysvisutilpy")
        toplot = np.atleast_1d(toplot).flatten()
        nplot = toplot.size

        _kw = [{} for i in range(nplot)]  # make ``nplot`` empty dicts.
        if plot_kw is not None:
            for k, vs in plot_kw.items():
                # If only one is given, duplicate for all ``toplot``:
                vs = np.atleast_1d(vs)
                if vs.size == 1:
                    vs = vs.tolist() * nplot
                for i, v in enumerate(vs):
                    _kw[i][k] = v

        for i, c in enumerate(toplot):
            ax.plot(self.datetime.jd, self.queried[c], **_kw[i])

        tax = yvu.append_xdate(ax, self.plotdate, self.queried[c], **kwargs)[0]
        return tax

    def highlight_time(self, ax, start, stop, format='isot', **kw):
        start = Time(start, format=format).jd
        stop = Time(stop, format=format).jd
        ax.fill_betweenx(ax.get_ylim(), start, stop, **kw)
        return (start, stop)

    # def anglesplot(
    #     self, ax,
    #     elmax_kw=dict(color='b', ls=':', r"Night $\mathrm{EL}_\mathrm{max}$"),
    #     elmax_window_day=6,
    #     alpha_kw=dict(color='b', ls='-', label=r"$\alpha$"),
    #     lunar_kw=dict(color='b', marker='o', alpha=0.5,
    #                   label='Lunar elongation\n(size=illumination)'),
    #     lunar_interval_day=6, lunar_marker_size_factor=1,
    #     ylabel_kw=dict(ylabel=r'Angles [$^\circ$]', color='b'),
    #     yticks_kw=dict(color='b')
    # ):
    #     try:
    #         import bottleneck as bn
    #     except ImportError:
    #         raise ImportError("Bottleneck is required for this."
    #                           + "https://github.com/pydata/bottleneck")

    #     eph = self.ephemerides()
    #     self.datetime = Time(eph['datetime_jd'], format='jd')
    #     self.elmax = bn.move_max(eph['EL'], window=24*elmax_window_day)
    #     plotdate = self.datetime.plot_date

    #     # ax.plot_date(np.nan, np.nan, color=elmax_kw['color'], , label='V-mag')
    #     lunar_sl = slice(None, None, 24*lunar_interval_day)
    #     ax.plot_date(plotdate, eph['alpha'], **alpha_kw)
    #     ax.plot_date(plotdate, elmax, **elmax_kw)
    #     ax.scatter(
    #         plotdate[lunar_sl],
    #         eph['lunar_elong'][lunar_sl],
    #         s=lunar_marker_size_factor*eph['lunar_illum'][lunar_sl],
    #         **lunar_kw
    #     )
    #     ax.set_ylabel(**ylabel_kw)
    #     ax.tick_params('y', **yticks_kw)
    #     # below are too easy to modify by the user,
    #     # so I didn't put freedom here.
    #     ax.set(
    #         xlabel='Date (UT)',
    #         ylim=(0, 180),
    #         title=self.shortname
    #     )
    #     ax.tick_params('y', colors=c1)
    #     ax.grid(axis='x')
    #     ax.legend(bbox_to_anchor=(1, -0.04), ncol=4)
    #     plt.gcf().autofmt_xdate()

    # ax2 = ax.twinx()
    # ax2.plot_date(datetime.plot_date[sl*7], eph['V'][sl*7], f'{c2}-')
    # ax2.set_ylabel(r'V-magnitude', color=c2)
    # ax2.tick_params('y', colors=c2)
