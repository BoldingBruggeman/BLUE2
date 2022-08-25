#!/usr/bin/env python

# Added the option to output at constant depths.


# was running with this:
#mpiexec -np $NPROCS --output-filename mpilog python -u $HOME/BLUE2/swes/swes_2022_06_22.py '2005-01-01 00:00:00' '2005-08-31 00:00:00' '/ACQUA/COMMONDATA/BLUE2-data/ERA5-europe/' --input_dir /ACQUA/COMMONDATA/BLUE2-data/swes/Input --tiling subdiv_96.pickle  --output_dir /BGFS/ACQUA/ferrnun/BLUE2/swes_gvc_2022_06_22/ --debug_output --save_restart /BGFS/ACQUA/ferrnun/BLUE2/swes_gvc_2022_06_22/restart.nc --initial



import argparse
import sys
import datetime
import os.path

import pygetm
import pygetm.legacy
import pygetm.input.tpxo


if os.path.isfile("../shared/blue2.py"):
    sys.path.insert(0, "../shared")
    import blue2
else:
    print("edit swes.py and set the folder where blue2.py is")
    quit()

setup = "swes"

parser = argparse.ArgumentParser()
blue2.config(parser)
args = parser.parse_args()

simstart = datetime.datetime.strptime(args.start, "%Y-%m-%d %H:%M:%S")
simstop = datetime.datetime.strptime(args.stop, "%Y-%m-%d %H:%M:%S")
if args.setup_dir is None:
    args.setup_dir = "."
if args.input_dir is None:
    args.input_dir = os.path.join(args.setup_dir, "Input")
if args.output is None:
    args.output_dir = os.path.join(args.setup_dir, ".")
tiling = args.tiling if args.tiling is not None else None
if args.meteo_dir is None:
    args.no_meteo = True
profile = "swes" if args.profile is not None else None

# Setup domain and simulation
domain = pygetm.legacy.domain_from_topo(
    os.path.join(args.setup_dir, "topo.nc"),
    nlev=45,
    z0_const=0.005,
    z0=0.005,
    Dmin=0.7,
    Dcrit=2.1,
    vertical_coordinate_method=pygetm.domain.VerticalCoordinates.GVC,
    ddu=2.0,
    ddl=1.0,
    Dgamma=50.0,
    Am=1.8e-06,
    # gamma_surf = True,
    tiling=args.tiling,
)

if args.boundaries:
    pygetm.legacy.load_bdyinfo(domain, os.path.join(args.setup_dir, "bdyinfo.dat"))
if args.rivers:
    pygetm.legacy.load_riverinfo(domain, os.path.join(args.setup_dir, "riverinfo.dat"))

if args.no_meteo:
    airsea = pygetm.airsea.Fluxes()
else:
    airsea=pygetm.airsea.FluxesFromMeteo(humidity_measure=pygetm.airsea.HumidityMeasure.DEW_POINT_TEMPERATURE, calculate_evaporation=True, longwave_method=pygetm.airsea.LongwaveMethod.JOSEY2, albedo_method = pygetm.airsea.AlbedoMethod.PAYNE)
    #airsea = (
    #    pygetm.airsea.FluxesFromMeteo(
    #        humidity_measure=pygetm.airsea.HumidityMeasure.DEW_POINT_TEMPERATURE
    #    ),
    #)

# Setup simulation
sim = pygetm.Simulation(
    domain,
    runtype=pygetm.BAROCLINIC,
    advection_scheme=pygetm.AdvectionScheme.UPSTREAM, #HSIMT,
    internal_pressure_method=pygetm.InternalPressure.SHCHEPETKIN_MCWILLIAMS,
#    internal_pressure_method=pygetm.InternalPressure.SHCHEPETKIN_MCWILLIAMS,
    airsea=airsea,
    #airsea=pygetm.airsea.FluxesFromMeteo(humidity_measure=pygetm.airsea.HumidityMeasure.DEW_POINT_TEMPERATURE),
    gotm=os.path.join(args.setup_dir, "gotmturb.nml"),
)

# sim.input_manager.debug_nc_reads()

# Open boundary conditions
if domain.open_boundaries:
    if args.tpxo9_dir is None:
        sim.logger.info("Reading 2D boundary data from file")
        domain.open_boundaries.z.set(pygetm.input.from_nc(os.path.join(args.input_dir, 'bdy_2d.nc'), 'elev'))
        domain.open_boundaries.u.set(pygetm.input.from_nc(os.path.join(args.input_dir, 'bdy_2d.nc'), 'u'))
        domain.open_boundaries.v.set(pygetm.input.from_nc(os.path.join(args.input_dir, 'bdy_2d.nc'), 'v'))
    else:
        sim.logger.info("Setting up TPXO tidal boundary forcing")
        TPXO9_data_dirs = ("/ACQUA/COMMONDATA")
        tpxo_dir = os.path.join("/ACQUA/COMMONDATA", "TPXO9")
        bdy_lon = domain.T.lon.all_values[domain.bdy_j, domain.bdy_i]
        bdy_lat = domain.T.lat.all_values[domain.bdy_j, domain.bdy_i]
        if domain.open_boundaries:
            sim.zbdy.set(
                pygetm.input.tpxo.get(bdy_lon, bdy_lat, root=tpxo_dir), on_grid=True
            )
            sim.bdyu.set(
                pygetm.input.tpxo.get(bdy_lon, bdy_lat, variable="u", root=tpxo_dir),
                on_grid=True,
            )
            sim.bdyv.set(
                pygetm.input.tpxo.get(bdy_lon, bdy_lat, variable="v", root=tpxo_dir),
                on_grid=True,
            )

    # Need to either support old GETM full field bdy handling - or extract only those points needed
    if True and sim.runtype == pygetm.BAROCLINIC:
        sim.temp.open_boundaries.type = pygetm.SPONGE
        sim.temp.open_boundaries.values.set(
            pygetm.input.from_nc(os.path.join(args.input_dir, "bdy_3d.nc"), "temp"),
            on_grid=True,
        )
        sim.salt.open_boundaries.type = pygetm.SPONGE
        sim.salt.open_boundaries.values.set(
            pygetm.input.from_nc(os.path.join(args.input_dir, "bdy_3d.nc"), "salt"),
            on_grid=True,
        )

# Initial salinity and temperature
if args.initial and sim.runtype == pygetm.BAROCLINIC:
    sim.logger.info("Setting up initial salinity and temperature conditions")
    sim.temp.set(
        pygetm.input.from_nc(
            os.path.join(args.input_dir, "initial.nc"), "temp"
        ),
        on_grid=True,
    )
    sim.salt.set(
        pygetm.input.from_nc(
            os.path.join(args.input_dir, "initial.nc"), "salt"
        ),
        on_grid=True,
    )
    sim.temp[..., domain.T.mask == 0] = pygetm.constants.FILL_VALUE
    sim.salt[..., domain.T.mask == 0] = pygetm.constants.FILL_VALUE
    sim.density.convert_ts(sim.salt, sim.temp)

if domain.rivers:
    for name, river in domain.rivers.items():
        river.flow.set(
            pygetm.input.from_nc(os.path.join(args.input_dir, "mergerivers4.nc"), name)
        )
        if sim.runtype == pygetm.BAROCLINIC:
#            river["salt"].set(0.01)
            pygetm.input.from_nc(
                os.path.join(args.input_dir, "mergerivers1.nc"), "%s_salt" % name
            )

#    import netCDF4
#    riverfile = os.path.join(args.input_dir, "mergerivers1.nc")
#    nc = netCDF4.Dataset(riverfile)
#
#    for name, river in domain.rivers.items():
#       #We can have several river-mouths x river:
#       river.flow.set(pygetm.input.from_nc(riverfile, river.original_name) /
#       river.split)
#       river["salt"].follow_target_cell = "%s_salt" % river.original_name not in nc.variables
#       if not river["salt"].follow_target_cell:
#          river["salt"].set(pygetm.input.from_nc(riverfile, "%s_salt" % river.original_name,))
#
#       river["temp"].follow_target_cell = "%s_temp" % river.original_name not in nc.variables
#       if not river["temp"].follow_target_cell:
#          river["temp"].set(pygetm.input.from_nc(riverfile, "%s_temp" % river.original_name,))

            

# Meteorological forcing - select between ERA-interim or ERA5
if not args.no_meteo:
    if False:
        sim.logger.info("Setting up ERA interim meteorological forcing")
        ERA_path = os.path.join(args.meteo_dir, "analysis_2015_%04i.nc" % simstart.year)
        era_kwargs = {"preprocess": lambda ds: ds.isel(time=slice(4, -4))}
    else:
        sim.logger.info("Setting up ERA5 meteorological forcing")
        ERA_path = os.path.join(args.meteo_dir, "era5_%04i.nc" % simstart.year)
        era_kwargs = {}
    sim.airsea.u10.set(pygetm.input.from_nc(ERA_path, "u10", **era_kwargs))
    sim.airsea.v10.set(pygetm.input.from_nc(ERA_path, "v10", **era_kwargs))
    sim.airsea.t2m.set(pygetm.input.from_nc(ERA_path, "t2m", **era_kwargs) - 273.15)
    sim.airsea.d2m.set(pygetm.input.from_nc(ERA_path, "d2m", **era_kwargs) - 273.15)
    sim.airsea.tcc.set(pygetm.input.from_nc(ERA_path, "tcc", **era_kwargs))
    sim.airsea.sp.set(pygetm.input.from_nc(ERA_path, "sp", **era_kwargs))
    sim.airsea.tp.set(pygetm.input.from_nc(ERA_path, 'tp', **era_kwargs) / 3600) 

if pygetm.BAROCLINIC:
    sim.radiation.set_jerlov_type(pygetm.radiation.JERLOV_II)

if sim.fabm:
    sim.logger.info("Setting up FABM dependencies that GETM does not provide")
    # sim.logger.info('Setting up FABM dependencies that GETM does not provide')
    # sim.fabm.get_dependency('downwelling_photosynthetic_radiative_flux').set(0)
    # sim.fabm.get_dependency('surface_downwelling_photosynthetic_radiative_flux').set(0)
    # sim.fabm.get_dependency('bottom_stress').set(0)
    # sim.fabm.get_dependency('bottom_stress').set(0)
    # if sim.runtype == pygetm.BAROTROPIC_3D:
    # sim.fabm.get_dependency('temperature').set(5.)
    # sim.fabm.get_dependency('practical_salinity').set(35.)

if args.output:
    sim.logger.info("Setting up output")


    if sim.runtype > pygetm.BAROTROPIC_2D:
        # -------------- File only averaging on MONTHS or YEARS
        #output = sim.output_manager.add_netcdf_file(
        #    os.path.join(args.output_dir, "medsea_monthly.nc"),
        #    interval=1,
        #    interval_units=pygetm.TimeUnit.MONTHS, #TimeUnit.YEARS
        #    sync_interval=True,
        #)
        #output.request(('temp', 'salt',), mask=True, time_average=True)
        #output.request(('uk', 'vk', 'ww',), grid=domain.T, z=pygetm.CENTERS, mask=True, time_average=True)
        # File with only temperature and salinity - daily
        # --------------- daily outputs
        output = sim.output_manager.add_netcdf_file(
            os.path.join(args.output_dir, "swes_temp_salt.nc"),
            interval=datetime.timedelta(hours=24),
            sync_interval=True,
        )
        output.request(("Ht",), mask=True)
    if sim.runtype == pygetm.BAROCLINIC:
        output.request(
            (
                "temp",
                "salt",
                "rho",
                "sst",
                "hnt",
                "qe",
                "qh", 
                "ql", 
                "swr", 
                "albedo", 
                "zen"
            ), 
            mask=True, 
            time_average=True
        )
        output.request(
            (
                "NN",
                "nuh",
                "SS",
                "num",
                "rad",
            ),
        )
        output = sim.output_manager.add_netcdf_file(
            os.path.join(args.output_dir, "swes_temp_salt_z.nc"),
            interval=datetime.timedelta(hours=24),
            sync_interval=True,
        )

        output.request("temp", "salt","uk","vk","ww", z=[-1000, -500, -100, -50, -3])
        if sim.fabm:
            output.request(("par", "med_ergom_o2", "med_ergom_OFL", "med_ergom_dd"),time_average=True)

if args.save_restart:
    sim.output_manager.add_restart(args.save_restart)

if args.load_restart:
    simstart = sim.load_restart(args.load_restart)

sim.start(
    simstart,
    timestep=13.0,
    split_factor=4,
    report=datetime.timedelta(hours=1),
    report_totals=datetime.timedelta(days=1),
    profile=profile,
)
while sim.time < simstop:
    sim.advance()
sim.finish()
