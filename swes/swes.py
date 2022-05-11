#!/usr/bin/env python

import argparse
import sys
import datetime
import os.path

import pygetm
import pygetm.legacy
import pygetm.input.tpxo

if os.path.isfile('../shared/blue2.py'):
    sys.path.insert(0, '../shared');
    import blue2
else:
    print('edit swes.py and set the folder where blue2.py is')
    quit()

setup = 'swes'

parser = argparse.ArgumentParser()
blue2.config(parser)
args = parser.parse_args()

simstart = datetime.datetime.strptime(args.start, '%Y-%m-%d %H:%M:%S')
simstop = datetime.datetime.strptime(args.stop, '%Y-%m-%d %H:%M:%S')
if args.setup_dir is None: args.setup_dir = '.'
if args.input_dir is None: args.input_dir = os.path.join(args.setup_dir,'input')
if args.output is None: args.output_dir = os.path.join(args.setup_dir,'.')
tiling = args.tiling if args.tiling is not None else None
if args.meteo_dir is None: args.no_meteo = True
profile = 'swes' if args.profile is not None else None

# Setup domain and simulation
domain = pygetm.legacy.domain_from_topo(os.path.join(args.setup_dir, 'topo.nc'), nlev=30, z0_const=0.001, tiling=args.tiling)
if args.boundaries:
    pygetm.legacy.load_bdyinfo(domain, os.path.join(args.setup_dir, 'bdyinfo.dat'))
if args.rivers:
    pygetm.legacy.load_riverinfo(domain, os.path.join(args.setup_dir, 'riverinfo.dat'))

if args.no_meteo:
    airsea=airsea=pygetm.airsea.Fluxes()
else:
    airsea=pygetm.airsea.FluxesFromMeteo(humidity_measure=pygetm.airsea.HumidityMeasure.DEW_POINT_TEMPERATURE),

# Setup simulation
sim = pygetm.Simulation(domain,
                        runtype=pygetm.BAROCLINIC,
#                        runtype=pygetm.BAROTROPIC,
                        advection_scheme=pygetm.AdvectionScheme.UPSTREAM,  # HSIMT
#                        airsea=airsea,
#                        airsea=pygetm.airsea.Fluxes(),
                        airsea=pygetm.airsea.FluxesFromMeteo(humidity_measure=pygetm.airsea.HumidityMeasure.DEW_POINT_TEMPERATURE),
#KB                        turbulence=pygetm.mixing.Turbulence(),
                        gotm=os.path.join(args.setup_dir, 'gotmturb.nml'),
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
if args.initial and sim.runtype == pygetm.BAROCLINIC:
    sim.logger.info('Setting up initial salinity and temperature conditions')
    sim.temp.set(pygetm.input.from_nc(os.path.join(args.input_dir, 'initial.nc'), 'temp'), on_grid=True) #KB
    sim.salt.set(pygetm.input.from_nc(os.path.join(args.input_dir, 'initial.nc'), 'salt'), on_grid=True) #KB
    sim.temp[..., domain.T.mask==0] = pygetm.constants.FILL_VALUE
    sim.salt[..., domain.T.mask==0] = pygetm.constants.FILL_VALUE
    sim.density.convert_ts(sim.salt, sim.temp)
else:
    sim.temp.set(15.)
    sim.salt.set(35.)

if domain.rivers:
    for name, river in domain.rivers.items():
        river.flow.set(pygetm.input.from_nc(os.path.join(args.input_dir, 'rivers.nc'), name))
        if sim.runtype == pygetm.BAROCLINIC:
            river['salt'].set(pygetm.input.from_nc(os.path.join(args.input_dir, 'rivers.nc'), '%s_salt' % name))

# Meteorological forcing - select between ERA-interim or ERA5
if not args.no_meteo:
    if False:
        sim.logger.info('Setting up ERA interim meteorological forcing')
        ERA_path = os.path.join(args.meteo_dir, 'analysis_2015_%04i.nc' % simstart.year)
        era_kwargs = {'preprocess': lambda ds: ds.isel(time=slice(4, -4))}
    else:
        sim.logger.info('Setting up ERA5 meteorological forcing')
        ERA_path = os.path.join(args.meteo_dir, 'era5_%04i.nc' % simstart.year)
        era_kwargs = {}
    sim.airsea.u10.set(pygetm.input.from_nc(ERA_path, 'u10', **era_kwargs))
    sim.airsea.v10.set(pygetm.input.from_nc(ERA_path, 'v10', **era_kwargs))
    sim.airsea.t2m.set(pygetm.input.from_nc(ERA_path, 't2m', **era_kwargs) - 273.15)
    sim.airsea.d2m.set(pygetm.input.from_nc(ERA_path, 'd2m', **era_kwargs) - 273.15)
    sim.airsea.tcc.set(pygetm.input.from_nc(ERA_path, 'tcc', **era_kwargs))
    sim.airsea.sp.set(pygetm.input.from_nc(ERA_path, 'sp', **era_kwargs))

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

if args.output:
    sim.logger.info('Setting up output')
    if not args.no_meteo:
        output = sim.output_manager.add_netcdf_file(os.path.join(args.output_dir, 'meteo.nc'), interval=datetime.timedelta(hours=1), sync_interval=None)
        output.request(('u10', 'v10', 'sp', 't2m', 'd2m', 'tcc'))
        if args.debug_output:
            output.request(('rhoa', 'qa', 'qs', 'qe', 'qh', 'ql', 'swr', 'albedo', 'zen'))
    output = sim.output_manager.add_netcdf_file(os.path.join(args.output_dir, 'swes_2d.nc'), interval=datetime.timedelta(hours=1), sync_interval=None)
    output.request(('Ht', ), mask=True)
    output.request(('zt', 'Dt', 'u1', 'v1', 'tausxu', 'tausyv', ))
    if args.debug_output:
        output.request(('U', 'V'), mask=True)
        output.request(('maskt', 'masku', 'maskv', ))   #, 'u_taus'
        output.request(('Du', 'Dv', 'dpdx', 'dpdy', 'z0bu', 'z0bv', 'z0bt'))   #, 'u_taus'
        output.request(('ru', 'rru', 'rv', 'rrv', ))
    if sim.runtype > pygetm.BAROTROPIC_2D:
        output = sim.output_manager.add_netcdf_file(os.path.join(args.output_dir, 'swes_3d.nc'), interval=datetime.timedelta(hours=6), sync_interval=None)
        output.request(('Ht', ), mask=True)
        output.request(('zt', 'uk', 'vk', 'ww', 'SS', 'num',))
        if args.debug_output:
            output.request(('fpk', 'fqk', 'advpk', 'advqk',))
    if sim.runtype == pygetm.BAROCLINIC:
        output.request(('temp', 'salt', 'rho', 'NN', 'rad', 'sst', 'hnt', 'nuh',))
        if args.debug_output:
            output.request(('idpdx', 'idpdy',))
        if sim.fabm_model:
            output.request(('par', 'med_ergom_o2', 'med_ergom_OFL', 'med_ergom_dd'))

sim.start(simstart, timestep=15., split_factor=20, report=180, profile=profile)
while sim.time < simstop:
    sim.advance()
sim.finish()
