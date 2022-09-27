#!/usr/bin/env python

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
    print("edit nwes.py and set the folder where blue2.py is")
    quit()

setup = "nwes"

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
profile = "nwes" if args.profile is not None else None

# Setup domain and simulation:
# https://pygetm.readthedocs.io/en/latest/api/pygetm.domain.html?highlight=domain
domain = pygetm.legacy.domain_from_topo(
    os.path.join(args.setup_dir, "topo.nc"),
    nlev=24,
    #z0=0.001,
    Dmin=5,
    Dcrit=10,
    vertical_coordinate_method=pygetm.domain.VerticalCoordinates.GVC,
    ddu=4.0,
    ddl=1.0,
    Dgamma=20.0,
    # gamma_surf = True,
    tiling=args.tiling,
)

if args.boundaries:
    pygetm.legacy.load_bdyinfo(domain, os.path.join(args.setup_dir, "bdyinfo.dat"))
if args.rivers:
    pygetm.legacy.load_riverinfo(domain, os.path.join(args.setup_dir, "GREENriver_19902018.dat"))

# Configure the airsea exchange:
# https://pygetm.readthedocs.io/en/latest/api/pygetm.airsea.html?highlight=airsea
if args.no_meteo:
    airsea = airsea = pygetm.airsea.Fluxes()
else:
    airsea = pygetm.airsea.FluxesFromMeteo(
        humidity_measure=pygetm.airsea.HumidityMeasure.DEW_POINT_TEMPERATURE,
        calculate_evaporation=True,
        longwave_method=pygetm.airsea.LongwaveMethod.JOSEY2,
    )

#if args.plot_domain:
#   f = domain.plot()
#   if f is not None:
#       f.savefig("domain_mesh.png")

# Configure the momentum class:
# https://pygetm.readthedocs.io/en/latest/api/pygetm.momentum.html?highlight=momentum
momentum = pygetm.momentum.Momentum(
    advection_scheme=pygetm.AdvectionScheme.SUPERBEE
)  # An=900., Am=0., pygetm.AdvectionScheme.SUPERBEE

# Setup simulation:
# https://pygetm.readthedocs.io/en/latest/api/pygetm.simulation.html?highlight=simulation
sim = pygetm.Simulation(
    domain,
    pygetm.BAROCLINIC,
    momentum=momentum,
    # Scalars
    advection_scheme=pygetm.AdvectionScheme.SUPERBEE,  # UPSTREAM SPLMAX13 SUPERBEE HSIMT P2_PDM MUSCL
    internal_pressure_method=pygetm.InternalPressure.SHCHEPETKIN_MCWILLIAMS,  # SHCHEPETKIN_MCWILLIAMS BLUMBERG_MELLOR
    airsea=airsea,
    gotm=os.path.join(args.setup_dir, "gotmturb.nml"),
    #fabm=os.path.join(args.setup_dir, "fabm.yaml"),
)
# Set spatial varying An - must be done after sim is defined
#sim.momentum.An.set(
#    pygetm.input.from_nc(os.path.join(args.input_dir, "m5a_An.nc"), "An"), on_grid=True
#)

#if args.plot_domain:
#    f = domain.plot(show_mesh=False, show_mask=True)
#    if f is not None:
#        f.savefig("domain_mask.png")

# sim.input_manager.debug_nc_reads()

# Open boundary conditions
if domain.open_boundaries:
    if args.tpxo9_dir is None:
        sim.logger.info("Reading 2D boundary data from file")
        domain.open_boundaries.z.set(pygetm.input.from_nc(os.path.join(args.input_dir, 'Boundary/bdy.2d.NWES.TPXO9OSU.1990.nc'), 'elev'))
        domain.open_boundaries.u.set(pygetm.input.from_nc(os.path.join(args.input_dir, 'Boundary/bdy.2d.NWES.TPXO9OSU.1990.nc'), 'u'))
        domain.open_boundaries.v.set(pygetm.input.from_nc(os.path.join(args.input_dir, 'Boundary/bdy.2d.NWES.TPXO9OSU.1990.nc'), 'v'))

# I think it should read the boundary as in hte old version, so I used that  


    else:
#        sim.logger.info("Setting up TPXO tidal boundary forcing")
#        TPXO9_data_dir = "/ACQUA/COMMONDATA"
#        tpxo_dir = os.path.join(TPXO9_data_dir, "TPXO9")
#        bdy_lon = domain.T.lon.all_values[domain.bdy_j, domain.bdy_i]
#        bdy_lat = domain.T.lat.all_values[domain.bdy_j, domain.bdy_i]
#        if domain.open_boundaries:
#            sim.zbdy.set(
#                pygetm.input.tpxo.get(bdy_lon, bdy_lat, root=tpxo_dir), on_grid=True
#            )
#            sim.bdyu.set(
#                pygetm.input.tpxo.get(bdy_lon, bdy_lat, variable="u", root=tpxo_dir),
#                on_grid=True,
#            )
#            sim.bdyv.set(
#                pygetm.input.tpxo.get(bdy_lon, bdy_lat, variable="v", root=tpxo_dir),
#                on_grid=True,
#            )
#  NOT WORKING : it doesnt like domain.bdy_j and stop at initialisation phase.. Not sure of which kind of data should be read (original TPXO9 ATlas ?) 




#    if sim.runtype == pygetm.BAROCLINIC:
#        sim.logger.info("Setting up temperature & salinity boundaries")
#        sim.temp.open_boundaries.type = pygetm.SPONGE  # CLAMPED
#        sim.salt.open_boundaries.type = pygetm.SPONGE  # CLAMPED
#        salt_bdy, temp_bdy = sim.density.lazy_convert_ts(
#            pygetm.input.from_nc(
#                os.path.join(args.input_dir, "Boundary/bdy.3d.NWA.RF_JvdM.v01.1990.nc"), "salt"
#            ),
#            pygetm.input.from_nc(
#                os.path.join(args.input_dir, "Boundary/bdy.3d.NWA.RF_JvdM.v01.1990.nc"), "temp"
#            ),
#            lon=domain.open_boundaries.lon.allgather()[:, None],
#            lat=domain.open_boundaries.lat.allgather()[:, None],
#            in_situ=True,
#        )
#        sim.temp.open_boundaries.values.set(
#            temp_bdy,
#            climatology=False,
#        )
#        sim.salt.open_boundaries.values.set(
#            salt_bdy,
#            climatology=False,
#        )
#   THIS PART IS NOT WORKING : It reads the file but MPI crash !! 



#
# Initial salinity and temperature
if args.initial and sim.runtype == pygetm.BAROCLINIC:
    sim.logger.info("Setting up initial salinity and temperature conditions")
    sim.temp.set(
        pygetm.input.from_nc(
            os.path.join(args.input_dir, "initial.physics.WOA.nc"),
            "temp",
        ),
        on_grid=True,
    )
    sim.salt.set(
        pygetm.input.from_nc(
            os.path.join(args.input_dir, "initial.physics.WOA.nc"),

            "salt",
        ),
        on_grid=True,
    )
    sim.temp[..., domain.T.mask == 0] = pygetm.constants.FILL_VALUE
    sim.salt[..., domain.T.mask == 0] = pygetm.constants.FILL_VALUE
    sim.density.convert_ts(sim.salt, sim.temp, in_situ=True)

# River handling:
split_rivers = True
if domain.rivers:
    if split_rivers:
        import netCDF4

        riverfile = os.path.join(args.input_dir, "GREENriver_19902018.nc")
        nc = netCDF4.Dataset(riverfile)

        for name, river in domain.rivers.items():
            river.flow.set(
                pygetm.input.from_nc(riverfile, river.original_name) / river.split
            )
            river["salt"].follow_target_cell = (
                "%s_salt" % river.original_name not in nc.variables
            )
            if not river["salt"].follow_target_cell:
                river["salt"].set(
                    pygetm.input.from_nc(
                        riverfile,
                        "%s_salt" % river.original_name,
                    )
                )

            river["temp"].follow_target_cell = (
                "%s_temp" % river.original_name not in nc.variables
            )
            if not river["temp"].follow_target_cell:
                river["temp"].set(
                    pygetm.input.from_nc(
                        riverfile,
                        "%s_temp" % river.original_name,
                    )
                )
    else:
        for name, river in domain.rivers.items():
            river.flow.set(
                pygetm.input.from_nc(os.path.join(args.input_dir, "GREENriver_19902018.nc"), name)
            )
            if sim.runtype == pygetm.BAROCLINIC:
                river["salt"].set(
                    pygetm.input.from_nc(
                        os.path.join(args.input_dir, "GREENriver_19902018.nc"), "%s_salt" % name
                    )
                )

    if (
        sim.fabm
    ):  # ATENCIO: fabm llegeix 'jrc_med_ergom_*', pero 'rivers.nc' te 'jrc_medergom_*'
        river["N3_n"].follow_target_cell = (
            "%s_N3_n" % river.original_name not in nc.variables
        )
        if not river["N3_n"].follow_target_cell:
            river["N3_n"].set(
                pygetm.input.from_nc(
                    riverfile, "%s_N3_n" % river.original_name
                )
            )

            river["N1_p"].follow_target_cell = (
                "%s_N1_p" % river.original_name not in nc.variables
            )
            if not river["N1_p"].follow_target_cell:
                river["N1_p"].set(
                    pygetm.input.from_nc(
                        riverfile, "%s_N1_p" % river.original_name
                    )
                )


# Meteorological forcing - select between ERA-interim or ERA5
if not args.no_meteo:
    if False:
        sim.logger.info("Setting up ERA interim meteorological forcing")
        ERA_path = os.path.join(args.meteo_dir, "analysis_2015_%04i.nc" % simstart.year)
        era_kwargs = {"preprocess": lambda ds: ds.isel(time=slice(4, -4))}
    else:
        sim.logger.info("Setting up ERA5 meteorological forcing")
        ERA_path = os.path.join(args.meteo_dir, "era5_????.nc")
        era_kwargs = {}
    sim.airsea.u10.set(pygetm.input.from_nc(ERA_path, "u10", **era_kwargs))
    sim.airsea.v10.set(pygetm.input.from_nc(ERA_path, "v10", **era_kwargs))
    sim.airsea.t2m.set(pygetm.input.from_nc(ERA_path, "t2m", **era_kwargs) - 273.15)
    sim.airsea.d2m.set(pygetm.input.from_nc(ERA_path, "d2m", **era_kwargs) - 273.15)
    sim.airsea.tcc.set(pygetm.input.from_nc(ERA_path, "tcc", **era_kwargs))
    sim.airsea.sp.set(pygetm.input.from_nc(ERA_path, "sp", **era_kwargs))
    sim.airsea.tp.set(pygetm.input.from_nc(ERA_path, "tp", **era_kwargs) * 2.7778e-04)

if pygetm.BAROCLINIC:
    sim.radiation.set_jerlov_type(pygetm.radiation.JERLOV_II)

if sim.fabm:
    sim.logger.info("Setting up initial FABM BIOGEOCHEM fields")
    #    sim["jrc_med_ergom_pp"].set(
    #        pygetm.input.from_nc(os.path.join(args.input_dir, "medsea_bio_init.nc"), "jrc_med_ergom_pp"), on_grid=True)

    sim["N3_n"].set(
        pygetm.input.from_nc(
            os.path.join(args.input_dir, "initial.ERSEM.RF.v07.nc"), "N3_n"
        ),
        on_grid=True,
    )
    sim["N1_p"].set(
        pygetm.input.from_nc(
            os.path.join(args.input_dir, "initial.ERSEM.RF.v07.nc"), "N1_p"
        ),
        on_grid=True,
    )
    
    # SI que funciona!
    sim.logger.info(" ")
    sim.logger.info(" ==> Setting up boundary FABM BIOGEOCHEM fields")
    sim["N3_n"].open_boundaries.type = pygetm.SPONGE
    sim["N3_n"].open_boundaries.values.set(
        pygetm.input.from_nc(
            os.path.join(args.input_dir, "bdy.3d.bio.nc"), "N3_n"
        ),
        climatology=True,
    )

    sim["N1_p"].open_boundaries.type = pygetm.SPONGE
    sim["N1_p"].open_boundaries.values.set(
        pygetm.input.from_nc(
            os.path.join(args.input_dir, "bdy.3d.bio.nc"), "N1_p"
        ),
        climatology=True,
    )

    sim.logger.info(" ")

# Configure output:
# https://pygetm.readthedocs.io/en/latest/api/pygetm.output.html?highlight=output
if args.output:
# and not args.dryrun:
    sim.logger.info("Setting up output")
    if not args.no_meteo:
        output = sim.output_manager.add_netcdf_file(
            os.path.join(args.output_dir, "meteo.nc"),
            interval=datetime.timedelta(hours=24),
            sync_interval=None,
        )
        output.request(("u10", "v10", "sp", "t2m", "d2m", "tcc", "tp"))
        if args.debug_output:
            output.request(
                "rhoa", "qa", "qs", "qe", "qh", "ql", "swr", "albedo", "zen", "pe"
            )
    output = sim.output_manager.add_netcdf_file(
        os.path.join(args.output_dir, "nwes_2d.nc"),
        interval=datetime.timedelta(hours=120),
        sync_interval=None,
    )
    output.request(("Ht",), mask=True)
    #    output.request(("An",), mask=True)

    output.request("zt", "Dt", "u1", "v1", "tausxu", "tausyv", mask=True, time_average=True)
    if args.debug_output:
        output.request(("U", "V"), mask=True)
        output.request("maskt", "masku", "maskv")
        output.request("Du", "Dv", "dpdx", "dpdy", "z0bu", "z0bv", "z0bt")
        output.request("ru", "rru", "rv", "rrv")

    if sim.runtype > pygetm.BAROTROPIC_2D:
        output = sim.output_manager.add_netcdf_file(
            os.path.join(args.output_dir, "nwes_3d.nc"),
            interval=datetime.timedelta(hours=120),
            sync_interval=None,
        )
        output.request(("Ht",), mask=True)
        output.request("zt", "uk", "vk", "ww", "SS", "num")
        if args.debug_output:
            output.request("fpk", "fqk", "advpk", "advqk")
            output.request("SxA", "SyA", "SxD", "SyD", "SxF", "SyF")
    if sim.runtype == pygetm.BAROCLINIC:
        output.request("temp", "salt", "rho", "NN", "rad", "sst", "hnt", "nuh", mask=True, time_average=True)
        if args.debug_output:
            output.request("idpdx", "idpdy", "SxB", "SyB")
        if sim.fabm:
            output.request("N3_n", "N1_p", "O2_O", "P1_Chl", "P2_Chl", "P3_Chl","P4_Chl", mask=True, time_average=True)

if args.save_restart:
# and not args.dryrun:
    sim.output_manager.add_restart(
        args.save_restart, interval=datetime.timedelta(days=7)
    )
#    sim.output_manager.add_restart(args.save_restart)

if args.load_restart: 
#and not args.dryrun:
    simstart = sim.load_restart(args.load_restart)

sim.start(
    simstart,
    timestep=12.0,
    split_factor=40,
    report=datetime.timedelta(hours=1),
    report_totals=datetime.timedelta(days=1),
    profile=profile,
)
#if args.dryrun:
#    sim.logger.info("Making a dryrun - skipping sim.advance()")
#else:
while sim.time < simstop:
    sim.advance()
sim.finish()
