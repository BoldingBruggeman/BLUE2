#!/usr/bin/env python

import argparse
import pathlib
import datetime
import os.path

import pygetm
import pygetm.legacy
import pygetm.input.tpxo

setup = 'medsea'

parser = argparse.ArgumentParser()
parser.add_argument('--setup_dir', type=pathlib.Path, help='Path to medsea configuration files', default='.')
parser.add_argument('--input_dir', type=pathlib.Path, help='Path to medsea input files', default='input' )
parser.add_argument('--meteo_dir', type=pathlib.Path, help='Path to meteo forcing files')
parser.add_argument('--tpxo9_dir', type=pathlib.Path, help='Path to TPXO9 configuration files')
parser.add_argument('--initial', action='store_true', help='Read initial salinity and temerature conditions from file')
parser.add_argument('--tiling', type=argparse.FileType('r'), help='Path to tiling pickle file')
parser.add_argument('--no_meteo', action='store_true', help='No meteo forcing')
parser.add_argument('--no_boundaries', action='store_true', help='No open boundaries')
parser.add_argument('--no_rivers', action='store_true', help='No river input')
parser.add_argument('--no_output', action='store_false', dest='output', help='Do not save any results to NetCDF')
parser.add_argument('--debug_output', action='store_true', dest='debug_output', help='Do not save any results to NetCDF')
args = parser.parse_args()


if args.setup_dir is None:
    args.setup_dir = '.'

if args.input_dir is None:
    args.input_dir = os.path.join(args.setup_dir,'input')

print(args)


# Setup domain and simulation
domain = pygetm.legacy.domain_from_topo(os.path.join(args.setup_dir, 'topo.nc'), nlev=30, z0_const=0.001, tiling=args.tiling)
if not args.no_boundaries:
    pygetm.legacy.load_bdyinfo(domain, os.path.join(args.setup_dir, 'bdyinfo.dat'))
if not args.no_rivers:
    pygetm.legacy.load_riverinfo(domain, os.path.join(args.setup_dir, 'riverinfo.dat'))

#airsea=pygetm.airsea.Fluxes()
sim = pygetm.Simulation(domain, 
                        runtype=pygetm.BAROCLINIC,
#                        runtype=pygetm.BAROTROPIC,
                        advection_scheme=pygetm.AdvectionScheme.HSIMT,
#                        airsea=pygetm.airsea.Fluxes(),
                        airsea=pygetm.airsea.FluxesFromMeteo(humidity_measure=pygetm.airsea.HumidityMeasure.DEW_POINT_TEMPERATURE),
#                        gotm=os.path.join(args.setup_dir, 'gotmturb.nml'),
                       )

#sim.input_manager.debug_nc_reads()

# Open boundary conditions
if domain.open_boundaries:
    if args.tpxo9_dir is None:
        sim.logger.info('Reading 2D boundary data from file')
#KB           domain.open_boundaries.z.set(pygetm.input.from_nc(os.path.join(args.input_dir, 'bdy_2d.nc'), 'elev'))
#KB           domain.open_boundaries.u.set(pygetm.input.from_nc(os.path.join(args.input_dir, 'Forcing/2D/bdy.2d.2006.nc'), 'u'))
#KB           domain.open_boundaries.v.set(pygetm.input.from_nc(os.path.join(args.input_dir, 'Forcing/2D/bdy.2d.2006.nc'), 'v'))
    else:
        sim.logger.info('Setting up TPXO tidal boundary forcing')
        TPXO9_data_dirs = ( '/server/data', '../../../igotm/data', '/ACQUA/COMMONDATA' )
        tpxo_dir = os.path.join(TPXO9_data_dir, 'TPXO9')
        bdy_lon = domain.T.lon.all_values[domain.bdy_j, domain.bdy_i]
        bdy_lat = domain.T.lat.all_values[domain.bdy_j, domain.bdy_i]
        if domain.open_boundaries:
            sim.zbdy.set(pygetm.input.tpxo.get(bdy_lon, bdy_lat, root=tpxo_dir), on_grid=True)
            sim.bdyu.set(pygetm.input.tpxo.get(bdy_lon, bdy_lat, variable='u', root=tpxo_dir), on_grid=True)
            sim.bdyv.set(pygetm.input.tpxo.get(bdy_lon, bdy_lat, variable='v', root=tpxo_dir), on_grid=True)

    # Need to either support old GETM full field bdy handling - or extract only those points needed
    if False and sim.runtype == pygetm.BAROCLINIC:
        sim.temp.open_boundaries.type = pygetm.SPONGE
        sim.temp.open_boundaries.values.set(pygetm.input.from_nc(os.path.join(args.input_dir, 'bdy_3d.nc'), 'temp'), on_grid=True)
        sim.salt.open_boundaries.type = pygetm.SPONGE
        sim.salt.open_boundaries.values.set(pygetm.input.from_nc(os.path.join(args.input_dir, 'bdy_3d.nc'), 'salt'), on_grid=True)

# Initial salinity and temperature
print(args.initial)
if args.initial and sim.runtype == pygetm.BAROCLINIC:
    sim.logger.info('Setting up initial salinity and temperature conditions')
    sim.temp.set(pygetm.input.from_nc(os.path.join(args.input_dir, 'medsea_5x5-clim-dec.nc'), 'temp'), on_grid=True)
    sim.salt.set(pygetm.input.from_nc(os.path.join(args.input_dir, 'medsea_5x5-clim-dec.nc'), 'salt'), on_grid=True)

if domain.rivers:
    for name, river in domain.rivers.items():
        river.flow.set(pygetm.input.from_nc(os.path.join(args.input_dir, 'rivers.nc'), name))
        if sim.runtype == pygetm.BAROCLINIC:
            river['salt'].set(pygetm.input.from_nc(os.path.join(args.input_dir, 'rivers.nc'), '%s_salt' % name))

# Meteorological forcing
if not args.no_meteo:
    ERA_data_dirs = ('/server/data', '../../../igotm/data', '/ACQUA/COMMONDATA/ECMWF')
    ERA_data_dir = '/server/data/ERA-interim-org'
    ERA_path = os.path.join(ERA_data_dir, 'analysis_2015.nc')
    print(ERA_path)
    sim.logger.info('Setting up ERA meteorological forcing')
    era_kwargs = {}
    era_kwargs = {'preprocess': lambda ds: ds.isel(time=slice(4, -4))}
    sim.airsea.tcc.set(pygetm.input.from_nc(ERA_path, 'tcc', **era_kwargs))
    sim.airsea.t2m.set(pygetm.input.from_nc(ERA_path, 't2m', **era_kwargs) - 273.15)
    #JRC sim.airsea.t2m.set(pygetm.input.from_nc(ERA_path, 't2', **era_kwargs) - 273.15)
    sim.airsea.d2m.set(pygetm.input.from_nc(ERA_path, 'd2m', **era_kwargs) - 273.15)
    #JRC sim.airsea.d2m.set(pygetm.input.from_nc(ERA_path, 'dev2', **era_kwargs) - 273.15)
    sim.airsea.sp.set(pygetm.input.from_nc(ERA_path, 'sp', **era_kwargs))
    #JRCsim.airsea.sp.set(pygetm.input.from_nc(ERA_path, 'slp', **era_kwargs))
    sim.airsea.u10.set(pygetm.input.from_nc(ERA_path, 'u10', **era_kwargs))
    sim.airsea.v10.set(pygetm.input.from_nc(ERA_path, 'v10', **era_kwargs))

if pygetm.BAROCLINIC:
    sim.radiation.set_jerlov_type(pygetm.radiation.JERLOV_II)

if sim.fabm_model:
    sim.logger.info('Setting up FABM dependencies that GETM does not provide')
    #sim.logger.info('Setting up FABM dependencies that GETM does not provide')
    #sim.get_fabm_dependency('downwelling_photosynthetic_radiative_flux').set(0)
    #sim.get_fabm_dependency('surface_downwelling_photosynthetic_radiative_flux').set(0)
    #sim.get_fabm_dependency('bottom_stress').set(0)
    #sim.get_fabm_dependency('bottom_stress').set(0)
    #if sim.runtype == pygetm.BAROTROPIC_3D:
        #sim.get_fabm_dependency('temperature').set(5.)
        #sim.get_fabm_dependency('practical_salinity').set(35.)

#if args.output:
    #sim.logger.info('Setting up output')
    #outputdir = '/ACQUA/BLUE22/blacksea'
    #outputdir = '/data/kb/BLUE2/out/blacksea'
    #outputdir = './out'
    #output = sim.output_manager.add_netcdf_file(os.path.join(outputdir, 'meteo.nc'), interval=1440) #1440
    #output.request(('u10', 'v10', 'sp', 't2m', 'd2m', 'tcc'))
    #output = sim.output_manager.add_netcdf_file(os.path.join(outputdir, '2d.nc'), interval=240) # 240
    #output.request(('U', 'V', 'zt'))
#
    #if runtype > 3:
        #output = sim.output_manager.add_netcdf_file(os.path.join(outputdir, '3d.nc'), interval=1440)
        #output.request(('dpdx', 'dpdy', 'tausxu', 'tausyv', 'u_taus'))
        #output.request(('temp', 'salt', 'nuh', 'NN'))
        #output.request(('med_ergom_o2', 'med_ergom_OFL'))

if args.output:
    sim.logger.info('Setting up output')
    #if not args.no_meteo:
        #output = sim.output_manager.add_netcdf_file('meteo.nc', interval=60, sync_interval=200000)
        #output.request(('u10', 'v10', 'sp', 't2m', 'qa', 'tcc'))
        #if args.debug_output:
            #output.request(('qe', 'qh', 'ql', 'swr', 'albedo', 'zen'))
    output = sim.output_manager.add_netcdf_file('medsea_2d.nc', interval=60, sync_interval=200000)
    output.request(('U', 'V'), mask=True)
    output.request(('zt', 'Dt', 'tausxu', 'tausyv', ))
    if args.debug_output:
        output.request(('maskt', 'masku', 'maskv', ))   #, 'u_taus'
        output.request(('Du', 'Dv', 'dpdx', 'dpdy', 'z0bu', 'z0bv', 'z0bt'))   #, 'u_taus'
        output.request(('ru', 'rru', 'rv', 'rrv'))
    if sim.runtype > pygetm.BAROTROPIC_2D:
        output = sim.output_manager.add_netcdf_file('medsea_3d.nc', interval=360, sync_interval=200000)
        output.request(('uk', 'vk', 'ww', 'SS', 'num',))
        if args.debug_output:
            output.request(('fpk', 'fqk', 'advpk', 'advqk',))
    if sim.runtype == pygetm.BAROCLINIC:
        output.request(('temp', 'salt', 'rho', 'NN', 'rad', 'sst', 'hnt', 'nuh',))
        if args.debug_output:
            output.request(( ))
        if sim.fabm_model:
            output.request(('par', 'med_ergom_o2', 'med_ergom_OFL', 'med_ergom_dd'))

sim.start(datetime.datetime(2015, 1, 1), timestep=10., split_factor=20, report=180)
#domain.update_depths() 
while sim.time < datetime.datetime(2015, 1, 1, 1):
    sim.advance()

sim.finish()


#BLUE2_setups_dir = next(filter(os.path.isdir, BLUE2_setups_dirs))
#TPXO9_data_dir = next(filter(os.path.isdir, TPXO9_data_dirs))
#ERA_data_dir = next(filter(os.path.isdir, ERA_data_dirs))
#ERA_path = os.path.join(ERA_data_dir, 'ERA-interim-org/analysis_2015.nc')
#JB ERA_path = os.path.join(ERA_data_dir, 'ERA-interim/2016.nc')
#ERA_path = os.path.join(ERA_data_dir, 'ncdf_erain/erain_2016_01.nc')

