import argparse
import pathlib

def config(parser):
#    parser = argparse.ArgumentParser()
    parser.add_argument('start', help='Similation start time - yyyy-mm-dd hh:mi:ss')
    parser.add_argument('stop', help='Similation stop time - yyyy-mm-dd hh:mi:ss')
    parser.add_argument('meteo_dir', type=pathlib.Path, help='Path to meteo forcing files')
    parser.add_argument('--setup_dir', type=pathlib.Path, help='Path to configuration files', default='.')
    parser.add_argument('--input_dir', type=pathlib.Path, help='Path to input files', default='Input' )
    parser.add_argument('--output_dir', type=pathlib.Path, help='Path to save output files', default='.')
    parser.add_argument('--tpxo9_dir', type=pathlib.Path, help='Path to TPXO9 configuration files')
    parser.add_argument('--tiling', help='Path to tiling pickle file')
    parser.add_argument('--initial', action='store_true', help='Read initial salinity and temerature conditions from file')
    parser.add_argument('--no_meteo', action='store_true', help='No meteo forcing')
    parser.add_argument('--no_boundaries', action='store_false', dest='boundaries', help='No open boundaries')
    parser.add_argument('--no_rivers', action='store_false', dest='rivers', help='No river input')
    parser.add_argument('--no_output', action='store_false', dest='output', help='Do not save any results to NetCDF')
    parser.add_argument('--debug_output', action='store_true', dest='debug_output', help='Do save additional fields to NetCDF')
    parser.add_argument('--profile', action='store_true', help='Save profiling data')
#    args = parser.parse_args()

#    if args.setup_dir is None: args.setup_dir = '.'
#    if args.input_dir is None: args.input_dir = os.path.join(args.setup_dir,'input')


