#!/usr/bin/env python3
"""Comparison TOPAZ - CTD within the Nordic Seas.

Python port of TOPAZ_NRT_CTD_Comp.m

Requirements:  pip install numpy scipy matplotlib netCDF4 gsw seawater

Usage:  python topaz_nrt_ctd_comp.py [YEAR MONTH]
        (prompts for year and month if they are not given)
"""

import sys
import warnings
from pathlib import Path

import gsw
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
import seawater                      # EOS-80, equivalent of sw_dens0
from scipy.spatial import cKDTree    # equivalent of dsearchn

# MATLAB's nanmean() silently returns NaN for empty/all-NaN input; numpy warns.
warnings.filterwarnings('ignore', category=RuntimeWarning)

# ----------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------
BASE = Path('/home/a23236/scratch/topas_tmp')   # NB: original had '/home//home/a23236/...'
COAST_FILE = BASE / 'Nordic4km.dat'
FIG_FORMATS = ('eps',)                          # export_fig '-eps'; add 'png', 'pdf' if wanted

MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
          'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
DYP = np.arange(5, 2001, 5)                     # interpolate to 5-m bins
DELTA_S = 0.03                                  # MLD density criterion

# Depth intervals (exclusive bounds, as in the original): 10-100, 100-300, 300-800, 800-2000 m
DD = [np.where((DYP > lo) & (DYP < hi))[0]
      for lo, hi in [(6, 101), (101, 301), (301, 801), (801, 2001)]]
DEPTH_LABELS = ['10-100 m', '100-300 m', '300-800 m', '800-2000 m']

RB_MAP = 'RdBu_r'                               # replaces redblue_colormap2 / rb_map
ASTERISK = (6, 2, 0)                            # MATLAB-style '*' marker

# Model sub-domain: MATLAB start [650 200], count [420 520] (1-based, x then y)
X0, NX, Y0, NY = 649, 420, 199, 520             # 0-based


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def nan_filled(a, dtype=float):
    """Masked array (netCDF fill values) -> ndarray with NaN, like MATLAB's ncread."""
    return np.ma.filled(np.ma.asarray(a).astype(dtype), np.nan)


def interp1(x, y, xi):
    """Linear interpolation like MATLAB interp1: NaN outside the data range."""
    ok = ~np.isnan(x)
    x, y = x[ok], y[ok]
    x, idx = np.unique(x, return_index=True)    # sorts and drops duplicate x
    return np.interp(xi, x, y[idx], left=np.nan, right=np.nan)


def dens0(s, t):
    """Density at surface pressure (sw_dens0)."""
    with np.errstate(invalid='ignore'):
        return seawater.dens0(s, t)


def conductivity_to_sm(c):
    """Original unit logic: <10 -> S/m, 10-50 -> mS/cm (divide by 10), negative or >50 -> NaN."""
    c = np.array(c, dtype=float)
    with np.errstate(invalid='ignore'):
        c[c < 0] = np.nan
        mid = (c > 10) & (c < 50)
        c[mid] = c[mid] / 10
        c[c > 50] = np.nan
    return c


def salinity(cond_sm, temp, depth):
    """Practical salinity from conductivity [S/m], temperature and depth/pressure.

    ASSUMPTION: stand-in for the original salinity.m (not provided). Uses TEOS-10/gsw
    (PSS-78) and treats depth as pressure in dbar. Check against your salinity.m.
    """
    return gsw.SP_from_C(cond_sm * 10.0, temp, depth)   # gsw wants mS/cm


def read_filelist(path):
    return Path(path).read_text().split()               # textscan(fid, '%s')


def depth_name(ds):
    return 'DEPH' if 'DEPH' in ds.variables else 'PRES'


def fillseg(ax, xy, face, edge):
    """Fill coastline polygons (segments separated by NaN rows)."""
    for seg in np.split(xy, np.flatnonzero(np.isnan(xy).any(axis=1))):
        seg = seg[~np.isnan(seg).any(axis=1)]
        if len(seg) > 2:
            ax.fill(seg[:, 0], seg[:, 1], facecolor=face, edgecolor=edge)


def save_fig(fig, stem):
    Path(stem).parent.mkdir(parents=True, exist_ok=True)
    for ext in FIG_FORMATS:
        fig.savefig(f'{stem}.{ext}', bbox_inches='tight', transparent=True)


# ----------------------------------------------------------------------------
# Load model data
# ----------------------------------------------------------------------------
def load_model(mfile, yr, mon):
    """Return lon, lat (x, y), depth, day number (days since 1970-01-01), T, S (x, y, z, t)."""
    with nc4.Dataset(mfile) as ds:
        # NetCDF files are (time, depth, y, x) in Python; MATLAB's ncread returns (x, y, depth, time),
        # so the 4-D fields are transposed to keep the original index order.
        # NB: assumes the file's dimension order is (time, depth, y, x).
        win = (slice(Y0, Y0 + NY), slice(X0, X0 + NX))
        m_lon = nan_filled(ds.variables['longitude'][win]).T
        m_lat = nan_filled(ds.variables['latitude'][win]).T
        m_dep = nan_filled(ds.variables['depth'][:])
        # Time is seconds since 1970-01-01 at 12:00 -> minus 0.5 d gives the calendar day
        m_day = np.rint(nan_filled(ds.variables['time'][:]) / 86400 - 0.5)
        m_tmp = nan_filled(ds.variables['thetao'][:, :, win[0], win[1]], np.float32).transpose(3, 2, 1, 0)
        m_sal = nan_filled(ds.variables['so'][:, :, win[0], win[1]], np.float32).transpose(3, 2, 1, 0)

    # Remove possible duplicates in time
    m_day, ta = np.unique(m_day, return_index=True)
    m_tmp = m_tmp[..., ta]
    m_sal = m_sal[..., ta]

    if yr == 2017 and mon == 11:     # Obvious error in field/day 19 in Nov 2017
        m_tmp[..., 18] = m_tmp[..., 17]
        m_sal[..., 18] = m_sal[..., 17]
    return m_lon, m_lat, m_dep, m_day, m_tmp, m_sal


# ----------------------------------------------------------------------------
# Load CTD data
# ----------------------------------------------------------------------------
def load_ctd(files):
    """Read all CTD files into one dict of arrays: (n_profiles,) or (n_profiles, n_levels)."""
    n_prof, n_lev = 0, 0
    for f in files:
        with nc4.Dataset(f) as ds:
            n_prof += len(ds.variables['TIME'])
            n_lev = max(n_lev, ds.variables[depth_name(ds)].shape[-1])

    obs = {k: np.full(n_prof, np.nan) for k in ('day', 'tqc', 'lon', 'lat', 'pqc')}
    obs.update({k: np.full((n_prof, n_lev), np.nan)
                for k in ('dep', 'depqc', 'temp', 'tempqc', 'sal', 'salqc')})

    k = 0
    for f in files:
        with nc4.Dataset(f) as ds:
            v = ds.variables
            n = len(v['TIME'])
            sl = slice(k, k + n)
            dname = depth_name(ds)                       # DEPH if present, otherwise PRES

            def prof(name):                              # (n_profiles, n_levels) block
                return np.atleast_2d(nan_filled(v[name][:]))

            dep = prof(dname)
            temp = prof('TEMP')
            dmp = dep.shape[1]

            # TIME is days since 1950-01-01; keep the calendar day (days since 1970-01-01)
            obs['day'][sl] = np.floor(nan_filled(v['TIME'][:])) - 7305
            obs['tqc'][sl] = nan_filled(v['TIME_QC'][:])
            obs['lon'][sl] = nan_filled(v['LONGITUDE'][:])
            obs['lat'][sl] = nan_filled(v['LATITUDE'][:])
            obs['pqc'][sl] = nan_filled(v['POSITION_QC'][:])
            obs['dep'][sl, :dmp] = dep
            obs['depqc'][sl, :dmp] = prof(dname + '_QC')
            obs['temp'][sl, :dmp] = temp
            obs['tempqc'][sl, :dmp] = prof('TEMP_QC')

            # Use CNDC if available, otherwise PSAL
            if 'CNDC' in v:
                cond = conductivity_to_sm(prof('CNDC'))
                obs['salqc'][sl, :dmp] = prof('CNDC_QC')
                obs['sal'][sl, :dmp] = salinity(cond, temp, dep)
            else:
                obs['sal'][sl, :dmp] = prof('PSAL')
                obs['salqc'][sl, :dmp] = prof('PSAL_QC')
            k += n

    # NB: as in the original, exact zeros are treated as missing (this also hits genuine
    # 0.0 values, e.g. 0.00 degC or QC flag 0).
    for a in obs.values():
        a[a == 0] = np.nan

    # Remove bad data (QC = 4)
    obs['temp'][obs['tempqc'] == 4] = np.nan
    obs['sal'][obs['salqc'] == 4] = np.nan

    # Remove unrealistic data that QC did not catch
    with np.errstate(invalid='ignore'):
        t_bad = (obs['temp'] > 40) | (obs['temp'] < -2)           # -2 < T < 40
        print('Ant. åpenbart feile temperaturobservasjoner tatt bort:')
        print(int(t_bad.sum()))
        obs['temp'][t_bad] = np.nan

        s_bad = (obs['sal'] > 50) | (obs['sal'] < 0)              # 0 < S < 50
        print('Ant. åpenbart feile saltholdighetsobservasjoner tatt bort:')
        print(int(s_bad.sum()))
        # NB: the original only counts these and does not remove them. To remove: obs['sal'][s_bad] = np.nan
    return obs


# ----------------------------------------------------------------------------
# Compare
# ----------------------------------------------------------------------------
def compare(obs, m_lon, m_lat, m_dep, m_day, m_tmp, m_sal):
    n_obs = obs['day'].size
    dc, dm = len(DYP), m_tmp.shape[2]

    def nan(*shape):
        return np.full(shape, np.nan)

    c_tmp_o, c_tmp_m = nan(n_obs, dc), nan(n_obs, dc)   # obs / model temp on common depths
    c_sal_o, c_sal_m = nan(n_obs, dc), nan(n_obs, dc)
    c_dns_o, c_dns_m = nan(n_obs, dc), nan(n_obs, dc)
    cm_tmp, cm_sal = nan(n_obs, dm), nan(n_obs, dm)     # model on model depths
    lon, lat = nan(n_obs), nan(n_obs)
    t_bias, s_bias = nan(4, n_obs), nan(4, n_obs)
    t_biast, s_biast, mld_bias = nan(n_obs), nan(n_obs), nan(n_obs)

    tree = cKDTree(np.column_stack([m_lon.ravel(), m_lat.ravel()]))
    n = 0                                               # 'teller'
    for t, day in enumerate(m_day):
        for i in np.where(obs['day'] == day)[0]:
            tqc, pqc = obs['tqc'][i], obs['pqc'][i]
            lon_o, lat_o = obs['lon'][i], obs['lat'][i]
            if (tqc < 1 or tqc > 2 or pqc < 1 or pqc > 2
                    or lon_o < -29 or lon_o > 49 or lat_o < 51 or lat_o > 81):
                continue
            if np.isnan(lon_o) or np.isnan(lat_o):      # (MATLAB would fail in dsearchn here)
                continue

            _, p_nh = tree.query([lon_o, lat_o])        # nearest model grid point
            im, jm = np.unravel_index(p_nh, m_lon.shape)

            ii = np.where(~np.isnan(obs['dep'][i]))[0]
            if ii.size == 0:
                continue

            lon[n], lat[n] = lon_o, lat_o
            dep = obs['dep'][i, ii]
            c_tmp_o[n] = interp1(dep, obs['temp'][i, ii], DYP)
            c_sal_o[n] = interp1(dep, obs['sal'][i, ii], DYP)
            c_dns_o[n] = dens0(c_sal_o[n], c_tmp_o[n])
            c_tmp_m[n] = interp1(m_dep, m_tmp[im, jm, :, t], DYP)
            c_sal_m[n] = interp1(m_dep, m_sal[im, jm, :, t], DYP)
            c_dns_m[n] = dens0(c_sal_m[n], c_tmp_m[n])
            cm_tmp[n] = m_tmp[im, jm, :, t]
            cm_sal[n] = m_sal[im, jm, :, t]

            # Bias within depth intervals (model - observation)
            for k, dk in enumerate(DD):
                t_bias[k, n] = np.nanmean(c_tmp_m[n, dk]) - np.nanmean(c_tmp_o[n, dk])
                s_bias[k, n] = np.nanmean(c_sal_m[n, dk]) - np.nanmean(c_sal_o[n, dk])

            # Bias over the whole water column 5-2000 m (levels where both have data)
            oo = np.where(~np.isnan(c_tmp_o[n]))[0]
            mm = np.where(~np.isnan(c_tmp_m[n]))[0]
            cc = np.intersect1d(oo, mm)
            t_biast[n] = np.nanmean(c_tmp_m[n, cc]) - np.nanmean(c_tmp_o[n, cc])
            s_biast[n] = np.nanmean(c_sal_m[n, cc]) - np.nanmean(c_sal_o[n, cc])

            # MLD bias
            with np.errstate(invalid='ignore'):
                ds1 = np.where(c_dns_o[n] - c_dns_o[n, 0] > DELTA_S)[0]
                ds2 = np.where(c_dns_m[n] - c_dns_m[n, 0] > DELTA_S)[0]
            if (ds1.size == 0 and ds2.size == 0) or cc.size == 0:
                mld_bias[n] = 0
            elif ds1.size == 0:                          # only model has a mixed-layer base
                mld_bias[n] = DYP[ds2.min()] - DYP[cc.max()]
            elif ds2.size == 0:                          # only obs has a mixed-layer base
                mld_bias[n] = DYP[cc.max()] - DYP[ds1.min()]
            else:
                mld_bias[n] = DYP[ds2.min()] - DYP[ds1.min()]
            n += 1

    res = dict(lon=lon, lat=lat, t_bias=t_bias, s_bias=s_bias, t_biast=t_biast,
               s_biast=s_biast, mld_bias=mld_bias, c_tmp_o=c_tmp_o, c_sal_o=c_sal_o,
               cm_tmp=cm_tmp, cm_sal=cm_sal)
    # Trim unused rows/columns: (4, n) bias arrays on the last axis, everything else on the first
    res = {k: a[:, :n] if k in ('t_bias', 's_bias') else a[:n] for k, a in res.items()}
    res['t_rmse'] = np.sqrt(np.nanmean(res['t_bias'] ** 2, axis=1))
    res['s_rmse'] = np.sqrt(np.nanmean(res['s_bias'] ** 2, axis=1))
    return res


# ----------------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------------
def plot_map(coast, lon, lat, values, vlim, title, xlabel_unit, figfile):
    fig, ax = plt.subplots()
    fillseg(ax, coast, '0.7', '0.7')
    sc = ax.scatter(lon, lat, s=36, c=values, cmap=RB_MAP, vmin=-vlim, vmax=vlim)
    cb = fig.colorbar(sc, ax=ax)
    if xlabel_unit:
        cb.ax.set_title(xlabel_unit, fontsize=12)
    ax.set_title(title, fontsize=12)
    ax.set_xlabel('Longitude [°E]', fontsize=12)
    ax.set_ylabel('Latitude [°N]', fontsize=12)
    ax.set_xlim(-25, 50)
    ax.set_ylim(54, 82)
    ax.tick_params(labelsize=12)
    save_fig(fig, figfile)


def plot_bias_depth(bias, rmse, xlim, xticks, text_x, xlabel, title, figfile):
    fig, ax = plt.subplots(figsize=(4, 6))
    for k in range(4):
        y = k + 1
        mean = np.nanmean(bias[k])
        ax.plot(bias[k], np.full(bias.shape[1], y), ls='none', marker=ASTERISK,
                ms=8, mew=1, color='0.7')
        ax.plot(mean, y, ls='none', marker=ASTERISK, ms=8, mew=1.5, color='k',
                label='Bias' if k == 0 else None)
        ax.plot(rmse[k], y, ls='none', marker=ASTERISK, ms=8, mew=1.5, color='r',
                label='RMSE' if k == 0 else None)
        ax.text(text_x[0], y + 0.15, f'{mean:.2f}', color='k', fontsize=12, va='center')
        ax.text(text_x[1], y + 0.15, f'{rmse[k]:.2f}', color='r', fontsize=12, va='center')
    ax.legend(ncol=2, loc='upper right')
    ax.set_title(title, fontsize=12)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_xlim(*xlim)
    ax.set_xticks(xticks)
    ax.set_ylim(4.5, 0.5)                               # ydir = reverse
    ax.set_yticks([1, 2, 3, 4])
    ax.set_yticklabels(DEPTH_LABELS)
    ax.tick_params(labelsize=12)
    save_fig(fig, figfile)


def plot_ts(res, title, figfile):
    ss = np.arange(32, 36.0001, 0.1)
    tt = np.arange(-2, 15.0001, 0.1)
    TT, SS = np.meshgrid(tt, ss)
    DD_ = dens0(SS, TT) - 1000

    fig, ax = plt.subplots()
    ax.plot(res['c_sal_o'].ravel(), res['c_tmp_o'].ravel(), '.', color='0.7', label='Obs')
    ax.plot(res['cm_sal'].ravel(), res['cm_tmp'].ravel(), 'r.', label='TOPAZ')
    ax.legend(loc='lower left')
    cs = ax.contour(SS, TT, DD_, levels=np.arange(20, 29.01, 0.5), colors='k', linestyles=':')
    ax.clabel(cs, fontsize=8)
    ax.set_title(title, fontsize=12)
    ax.set_xlabel('Salinity', fontsize=12)
    ax.set_ylabel('Temperature [deg C]', fontsize=12)
    ax.set_xlim(33, 35.5)
    ax.set_ylim(-2, 12)
    ax.tick_params(labelsize=12)
    save_fig(fig, figfile)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    if len(sys.argv) >= 3:
        yr, mon = int(sys.argv[1]), int(sys.argv[2])
    else:
        yr = int(input('Aar:'))
        mon = int(input('Mnd:'))

    fpath = BASE / f'figs{yr}' / f'M{mon:02d}'
    mpath = BASE / 'CTD_Valid' / str(yr) / f'M{mon}'
    mfile = mpath / f'TOPAZ_NRT_{yr}{mon:02d}.nc'
    flist = mpath / 'filelist.txt'

    m_lon, m_lat, m_dep, m_day, m_tmp, m_sal = load_model(mfile, yr, mon)
    obs = load_ctd(read_filelist(flist))
    res = compare(obs, m_lon, m_lat, m_dep, m_day, m_tmp, m_sal)

    coast = np.loadtxt(COAST_FILE)[:, :2]
    tag = f'{MONTHS[mon - 1]}{yr}'
    stem = lambda name: fpath / f'TOPAZ_CTD_Comparison_{name}_{yr}{mon:02d}'

    # Maps of bias
    plot_map(coast, res['lon'], res['lat'], res['t_biast'], 2,
             f'Temperature bias {tag}', '°C', stem('MapTemp'))
    plot_map(coast, res['lon'], res['lat'], res['s_biast'], 1,
             f'Salinity bias {tag}', None, stem('MapSalt'))
    plot_map(coast, res['lon'], res['lat'], res['mld_bias'], 250,
             f'MLD bias {tag}', 'm', stem('MapMLD'))

    # Bias vs. depth
    plot_bias_depth(res['t_bias'], res['t_rmse'], (-2, 2), np.arange(-2, 2.01, 1), (-1.5, 1),
                    'Temperature difference [deg C]', f'TOPAZ vs. CTD {tag}', stem('Temp'))
    plot_bias_depth(res['s_bias'], res['s_rmse'], (-1, 1), np.arange(-1, 1.01, 0.5), (-0.25, 0.25),
                    'Salinity difference', f'TOPAZ vs. CTD {tag}', stem('Salt'))

    # T-S diagram
    plot_ts(res, f'TOPAZ vs. CTD {tag}', stem('TS'))

    plt.show()


if __name__ == '__main__':
    main()