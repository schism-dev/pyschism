#! /usr/bin/env python
import argparse
import pathlib
from functools import lru_cache
from scipy.interpolate import griddata, interp1d
import matplotlib.pyplot as plt
from pyproj import Transformer
from datetime import datetime
import scipy.io as sio
import numpy as np
from datetime import timedelta
import imageio
import os
from pyschism.mesh import Hgrid


wt = os.getenv('TCPROF_WIND_TYPE')
wt = 'GAHM' if wt is None else wt


PARENT = pathlib.Path(__file__).parent.resolve()
FILE = PARENT / f'TC_Prof_all_weighted_composite_R2_{wt}-v6.mat'


class TCProf:

    def __init__(self, event_id):
        self._event_id = event_id

    def interpolate(
        self,
        hgrid,
        output_path,
        timestep,
        nskip=1,
        overwrite=False
    ):

        # certify hgrid
        msg = f"hgrid must be of type {Hgrid}"
        assert isinstance(hgrid, Hgrid), msg
        assert isinstance(nskip, int)

        xy = hgrid.get_xy(crs="EPSG:4326")

        def interp_field(step, x, y, field):
            return griddata(
                (x(step), y(step)),
                field(step),
                (xy[:, 0], xy[:, 1]),
                method='linear',
                fill_value=0.
                )

        # certify output_path
        output_path = pathlib.Path(output_path)
        if output_path.is_file() and not overwrite:
            msg = f"File {output_path} exists and overwrite=False"
            raise IOError(msg)

        steps, x, y, w, d, p = self._interp_hgrid(hgrid, timestep, nskip)

        with open(output_path, 'w') as f:
            for step in steps:
                line = f"{step.total_seconds()}, "
                wdata = interp_field(step, x, y, w)
                ddata = interp_field(step, x, y, d)
                pdata = interp_field(step, x, y, p)
                for i in range(hgrid.values.size):
                    u = wdata[i]*np.cos(ddata[i])
                    v = wdata[i]*np.sin(ddata[i])
                    line += f"{u} {v} {pdata[i]}, "
                line = line[:-2] + '\n'
                f.write(line)

    def plot_wfield(self, iframe, dsfact=6, scale=5000):
        fig, ax = plt.subplots(figsize=(16, 12))
        values = self.wind_magnitude[iframe]
        values = np.ma.masked_where(values < 64., values)
        ax_wind = ax.quiver(
            self.lon_arrays[iframe][::dsfact, ::dsfact],
            self.lat_arrays[iframe][::dsfact, ::dsfact],
            self.u_arrays[iframe][::dsfact, ::dsfact],
            self.v_arrays[iframe][::dsfact, ::dsfact],
            values[::dsfact, ::dsfact],
            scale=scale,
            cmap='jet'
            )
        ax.axis('scaled')
        fig.suptitle('Interpolated Wind Field of TC: ' + self.event_id)
        ax.set_title(self.datetime[iframe].strftime('%b %d, %Y %H:%M'))
        ax.set_xlabel('Longitude (°E)')
        ax.set_ylabel('Latitude (°N)')
        cbar = plt.colorbar(ax_wind, format='%.1f')
        cbar.ax.set_ylabel('Wind speed [km/h]', rotation=90)
        plt.show()
        plt.close(plt.gcf())

    def plot_pfield(self, iframe):
        fig, ax = plt.subplots(figsize=(16, 12))
        # ax.plot(self.coast[:, 1], self.coast[:, 0], color='black',
        #         linewidth=0.5)
        values = self.pressure[iframe]
        values = np.ma.masked_where(self.wind_magnitude[iframe] < 64., values)
        ax_press = ax.contourf(
            self.lon_arrays[iframe],
            self.lat_arrays[iframe],
            values,
            levels=256,
            cmap='jet',
            # vmax=self.plot_vmax,
            # vmin=self.plot_vmin
            )
        ax.axis('scaled')
        fig.suptitle('Interpolated Pressure Field of TC: ' + self.event_id)
        ax.set_title(self.datetime[iframe].strftime('%b %d, %Y %H:%M'))
        ax.set_xlabel('Longitude (°E)')
        ax.set_ylabel('Latitude (°N)')
        cbar = plt.colorbar(ax_press, format='%.1f')
        cbar.ax.set_ylabel('Pressure [hPa]', rotation=90)
        plt.show()
        plt.close(plt.gcf())

    def make_gif(self):
        images = list()
        dirname = self.filename
        os.makedirs(dirname, exist_ok=True)
        for i, time in enumerate(self.datetime):
            fig, ax = plt.subplots(figsize=(16, 12))
            ax.plot(self.coast[:, 1], self.coast[:, 0], color='black',
                    linewidth=0.5)
            values = self.values[i]
            values = np.ma.masked_where(values < 64., values)
            ax_wind = ax.contourf(self.lon_arrays[i], self.lat_arrays[i],
                                  values, levels=512, cmap='jet',
                                  vmax=self.plot_vmax, vmin=self.plot_vmin
                                  )
            ax.axis('scaled')
            fig.suptitle('Interpolated Wind Field of TC: ' + self.name[i])
            ax.set_title(time.strftime('%b %d, %Y %H:%M'))
            ax.set_xlabel('Longitude (°E)')
            ax.set_ylabel('Latitude (°N)')
            cbar = plt.colorbar(ax_wind, format='%.1f')
            cbar.ax.set_ylabel('Wind speed km/h', rotation=90)

            plt.savefig('{}/{:d}.png'.format(dirname, i))
            images.append(imageio.imread('{}/{:d}.png'.format(dirname, i)))
            plt.close()
        imageio.mimsave('{}/full_movie.gif'.format(dirname), images, fps=2)

    @property
    @lru_cache
    def tcprof(self):
        tcprof = {}
        for i, lines in enumerate(sio.loadmat(FILE)['TC_Prof'][0, :]):
            name = lines[0][0, 0][6][0]
            year = lines[0][0, 0][8][0]
            storm = (name + year).capitalize()
            tcprof[storm] = {
                "wind": [],
                "press": [],
                "param": [],
                "radio": [],
                "datetime": [],
            }
            for entry in lines[0][0, :]:
                if np.any(np.isnan(entry[3])):
                    continue
                tcprof[storm]['wind'].append(entry[0])
                tcprof[storm]['press'].append(entry[1])
                tcprof[storm]['param'].append(entry[3])
                tcprof[storm]['radio'].append(entry[5])
                tcprof[storm]['datetime'].append(
                    datetime.strptime(
                        f"{entry[9][0]} {entry[8][0]}",
                        "%b %d %H:%M %Y")
                    )
            tcprof[storm]['storm_id'] = f'AL{int(entry[7]):02d}{year}'
            tcprof[storm]['matlab_id'] = f'{int(i)+1}'
        return tcprof

    @property
    def event_id(self):
        return self._event_id

    @property
    def datetime(self):
        return self.tcprof[self.event_id]['datetime']

    @property
    @lru_cache
    def wind_magnitude(self):
        wind_magnitude = list()
        for date in self.datetime:
            wind_magnitude.append(self.wfields[date])
        return wind_magnitude

    @property
    @lru_cache
    def pressure(self):
        press = list()
        for date in self.datetime:
            press.append(self.pfields[date])
        return press

    @property
    @lru_cache
    def wfields(self):
        r_out, theta_out = np.meshgrid(
                np.linspace(0, 500000, 501),
                np.linspace(0., 2.*np.pi, 501)
                )
        fields = {}
        for i, date in enumerate(self.datetime):
            w_input = list()
            for j, angle in enumerate(np.linspace(-np.pi/4, 9*np.pi/4, num=6)):
                w_input.append(
                    self.tcprof[self.event_id]['wind'][i][((j-1) % 4), :])
            w_input = np.array(w_input).flatten()
            wfield = griddata(
                (self.r_input, self.theta_input),
                w_input,
                (r_out.flatten(), theta_out.flatten()),
                method='linear')
            # wfield interp
            idx = np.where(np.isnan(wfield))[0]
            _idx = np.where(~np.isnan(w_input))[0]
            _wfield = griddata(
                (self.r_input[_idx], self.theta_input[_idx]),
                w_input[_idx],
                (r_out.flatten()[idx], theta_out.flatten()[idx]),
                method='nearest')
            for i, _idx in enumerate(idx):
                wfield[_idx] = _wfield[i]
            wfield = wfield.reshape(r_out.shape)
            fields[date] = wfield
        return fields

    @property
    @lru_cache
    def pfields(self):
        r_out, theta_out = np.meshgrid(
                np.linspace(0, 500000, 501),
                np.linspace(0., 2.*np.pi, 501)
                )
        fields = {}
        for i, date in enumerate(self.datetime):
            p_input = list()
            for j, angle in enumerate(np.linspace(-np.pi/4, 9*np.pi/4, num=6)):
                p_input.append(
                        self.tcprof[self.event_id]['press'][i][((j-1) % 4), :],
                        )
            p_input = np.array(p_input).flatten()
            pfield = griddata(
                (self.r_input, self.theta_input),
                p_input,
                (r_out.flatten(), theta_out.flatten()),
                method='linear')
            # pfield interp
            idx = np.where(np.isnan(pfield))
            if len(idx[0]) > 0:
                _idx = np.where(~np.isnan(p_input))
                _pfield = griddata(
                    (self.r_input[_idx], self.theta_input[_idx]),
                    p_input[_idx],
                    (r_out.flatten()[idx], theta_out.flatten()[idx]),
                    method='nearest')
                for i, _idx in enumerate(idx):
                    pfield[_idx] = _pfield[i]
            pfield = pfield.reshape(r_out.shape)
            fields[date] = pfield
        return fields

    @property
    @lru_cache
    def wdirfields(self):
        r_out, theta_out = np.meshgrid(
                np.linspace(0, 500000, 501),
                np.linspace(0., 2.*np.pi, 501)
                )
        fields = {}
        for i, date in enumerate(self.datetime):
            wdir_input = list()
            for j, angle in enumerate(np.linspace(-np.pi/4, 9*np.pi/4, num=6)):
                wdir_input.append(
                    np.asarray(
                        501 * [np.pi / 2 + angle + self.wind_deviation_angle])
                    )
            wdir_input = np.array(wdir_input).flatten()
            wdirfield = griddata(
                (self.r_input, self.theta_input),
                wdir_input,
                (r_out.flatten(), theta_out.flatten()),
                method='linear')
            # wdirfield interp
            idx = np.where(np.isnan(wdirfield))
            if len(idx[0]) > 0:
                _idx = np.where(~np.isnan(wdir_input))
                _wdirfield = griddata(
                    (self.r_input[_idx], self.theta_input[_idx]),
                    wdir_input[_idx],
                    (r_out.flatten()[idx], theta_out.flatten()[idx]),
                    method='nearest')
                for i, _idx in enumerate(idx):
                    wdirfield[_idx] = _wdirfield[i]
            wdirfield = wdirfield.reshape(r_out.shape)
            fields[date] = wdirfield
        return fields

    @property
    @lru_cache
    def r_input(self):
        r_input = []
        for j, angle in enumerate(np.linspace(-np.pi/4, 9*np.pi/4, num=6)):
            r_input.append(np.linspace(0, 500000, 501))
        return np.array(r_input).flatten()

    @property
    @lru_cache
    def theta_input(self):
        theta_input = []
        for j, angle in enumerate(np.linspace(-np.pi/4, 9*np.pi/4, num=6)):
            theta_input.append(np.asarray(501 * [angle]))
        return np.array(theta_input).flatten()

    @property
    @lru_cache
    def lon_arrays(self):
        lon_arrays = []
        for date, field in self.lonlatfields.items():
            lon_arrays.append(field['lon'])
        return lon_arrays

    @property
    @lru_cache
    def lat_arrays(self):
        lat_arrays = []
        for date, field in self.lonlatfields.items():
            lat_arrays.append(field['lat'])
        return lat_arrays

    @property
    @lru_cache
    def u_arrays(self):
        u_arrays = []
        for date in self.datetime:
            u_arrays.append(self.wfields[date] * np.cos(self.wdirfields[date]))
        return u_arrays

    @property
    @lru_cache
    def v_arrays(self):
        v_arrays = []
        for date in self.datetime:
            v_arrays.append(self.wfields[date] * np.sin(self.wdirfields[date]))
        return v_arrays

    @property
    def wind_deviation_angle(self):
        try:
            return self.__wind_deviation_angle
        except AttributeError:
            return np.pi / 9.

    @property
    @lru_cache
    def wind_directions(self):
        wind_directions = list()
        for date in self.datetime:
            wind_directions.append(self.wdirfields[date])
        return wind_directions

    @property
    def eye_lon(self):
        try:
            return self.__eye_lon
        except AttributeError:
            self._init_eye_data()
            return self.__eye_lon

    @property
    def eye_lat(self):
        try:
            return self.__eye_lat
        except AttributeError:
            self._init_eye_data()
            return self.__eye_lat

    @property
    def eye_x(self):
        try:
            return self.__eye_x
        except AttributeError:
            self._init_eye_data()
            return self.__eye_x

    @property
    def eye_y(self):
        try:
            return self.__eye_y
        except AttributeError:
            self._init_eye_data()
            return self.__eye_y

    @property
    @lru_cache
    def lonlatfields(self):
        # initialize coordinate transformation objects
        Mercator_to_WGS84 = Transformer.from_proj(
            'EPSG:3395',
            'EPSG:4326',
            always_xy=True,
            )
        # precompute output base meshgrid in cartesian coordinates
        x, y = pol2cart(*np.meshgrid(
                np.linspace(0, 500000, 501),
                np.linspace(0, 2.*np.pi, 501)
                ))
        # add eye location to the base grid at each datetime
        fields = {}
        for i, date in enumerate(self.datetime):
            lon, lat = Mercator_to_WGS84.transform(
                self.eye_x[i] + x,
                self.eye_y[i] + y)
            fields[date] = {'lon': lon, 'lat': lat}
        return fields

    @wind_deviation_angle.setter
    def wind_deviation_angle(self, wind_deviation_angle):
        self.__wind_deviation_angle = float(wind_deviation_angle)

    def _init_eye_data(self):
        eye_lon = []
        eye_lat = []
        for i, date in enumerate(self.datetime):
            eye_lon.append(self.tcprof[self.event_id]['param'][i][9, 1])
            eye_lat.append(self.tcprof[self.event_id]['param'][i][9, 0])
        # precompute eye location in mercator coordinates
        WGS84_to_Mercator = Transformer.from_proj(
            'EPSG:4326',
            'EPSG:3395',
            always_xy=True,
            )
        eye_x, eye_y = WGS84_to_Mercator.transform(eye_lon, eye_lat)
        self.__eye_x = eye_x
        self.__eye_y = eye_y
        self.__eye_lon = eye_lon
        self.__eye_lat = eye_lat

    @property
    def _event_id(self):
        return self.__event_id

    @_event_id.setter
    def _event_id(self, event_id):
        msg = f'event_id must be one of {list(self.tcprof.keys())}'
        assert event_id.capitalize() in self.tcprof, msg
        self.__event_id = event_id

    def _interp_hgrid(self, hgrid, timestep, nskip=1):

        # get list of steps
        wtiminc = timedelta(seconds=timestep*nskip)
        total_time = self.datetime[-1] - self.datetime[0]
        current = timedelta(seconds=0)
        steps = []
        while current.total_seconds() < total_time.total_seconds():
            steps.append(current)
            current += wtiminc
        if steps[-1].total_seconds() < total_time.total_seconds():
            steps.append(current)

        # create "times" vector
        times = [(t - self.datetime[0]).total_seconds() for t in self.datetime]

        # create x, y interpolators
        x = interp1d(times, self.lon_arrays, axis=0)
        y = interp1d(times, self.lat_arrays, axis=0)

        # create pressure interpolator
        p = interp1d(times, self.pressure, axis=0)

        # create the w, d interpolator
        w = interp1d(times, self.wind_magnitude, axis=0)
        d = interp1d(times, self.wind_directions, axis=0)

        return steps, x, y, w, d, p


# auxilliary functions
def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)


# interface functions
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('event_id')
    parser.add_argument('--overwrite', action='store_true')

    # init subparsers
    subparsers = parser.add_subparsers(dest='mode')
    subparsers.required = True

    # data plotting subparsers
    plot_data = subparsers.add_parser('plot_data')
    plot_data.add_argument('plot_type', choices=['wind', 'press'])
    plot_data.add_argument('--index', type=int)

    # mesh plotting subparsers
    plot_mesh = subparsers.add_parser('plot_mesh')
    plot_mesh.add_argument('plot_type', choices=['wind', 'press'])
    plot_mesh.add_argument('--index', type=int)

    # anim subparsers
    save_anim = subparsers.add_parser('save_anim')
    save_anim.add_argument('plot_type', choices=['wind', 'press'])

    # interp subparsers
    interp = subparsers.add_parser('interp')
    interp.add_argument('timestep', type=float)
    interp.add_argument('hgrid')
    interp.add_argument('crs')
    interp.add_argument('output_path',)
    interp.add_argument('--nskip', type=int, default=1)

    # parse args
    args = parser.parse_args()

    return args


def interp(args, tcprof):
    if args.output_path is None:
        output_path = (pathlib.Path('.') / 'wind.th').resolve()
    else:
        output_path = pathlib.Path(args.output_path).resolve()
    tcprof.interpolate(
        Hgrid.open(args.hgrid, crs=args.crs),
        output_path,
        args.timestep,
        nskip=args.nskip,
        overwrite=args.overwrite
        )


def plot(args, tcprof):
    if args.plot_type == 'wind':
        if args.index is not None:
            tcprof.plot_wfield(args.index)
        else:
            msg = "Need to implement wind movie."
            raise NotImplementedError(msg)
    elif args.plot_type == 'press':
        if args.index is not None:
            tcprof.plot_pfield(args.index)
        else:
            msg = "Need to implement pressure movie."
            raise NotImplementedError(msg)


def save_anim(args, tcprof):
    if args.output_path is None:
        output_path = (
            pathlib.Path('.') / f'{args.anim_type}.gif').resolve()
    else:
        output_path = pathlib.Path(output_path).resolve()

    if args.anim_type == 'wind':
        tcprof.make_anim_wfield(args.output_directory)


def mesh_plot(args, tcprof):
    msg = 'Need to implement mesh plotting'
    raise NotImplementedError(msg)


def main():
    args = parse_args()
    tcprof = TCProf(args.event_id)

    # mesh interpolation mode
    if args.mode == 'interp':
        interp(args, tcprof)

    # tcprof data plot mode
    elif args.mode == 'plot_data':
        plot(args, tcprof)

    # tcprof gif generation mode
    elif args.mode == 'save_anim':
        save_anim(args, tcprof)

    # result visualization mode
    elif args.mode == 'plot_mesh':
        mesh_plot(args, tcprof)

    # duck-type error (should be unreachable)
    else:
        msg = f"Unrecognized mode: {args.mode}."
        raise Exception(msg)

    return 0


def init():
    if __name__ == '__main__':
        try:
            import colored_traceback
            colored_traceback.add_hook(always=True)
        except ModuleNotFoundError:
            pass
        exit(main())


init()
