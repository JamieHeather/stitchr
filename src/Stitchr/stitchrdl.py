# -*- coding: utf-8 -*-

"""
stitchrdl.py

Download IMGT/GENE-DB data for stitching, for a given species
"""


import os
from . import stitchrfunctions as fxn
import argparse
import subprocess
import sys
import shutil

# Make sure IMGTgeneDL is installed
try:
    import IMGTgeneDL
except (ImportError, ModuleNotFoundError) as err:
    print("IMGTgeneDL is required to run this script (pip install IMGTgeneDL)")
    sys.exit()


__version__ = '0.1.3'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'

sys.tracebacklimit = 0  # comment when debugging


def args():
    """
    args(): Obtains command line arguments which dictate the script's behaviour
    """

    # Help flag
    parser = argparse.ArgumentParser(
        description="stitchrdl v" + str(__version__) + '\n' + ": Process input data for use with stitchr. \n"
                    "See https://github.com/JamieHeather/stitchr and https://doi.org/10.1093/nar/gkac190.")

    # Input and output options
    parser.add_argument('-s', '--species', required=False, type=str, default='HUMAN',
                        help="Species (common name). Optional. Default = HUMAN")

    parser.add_argument('--version', action='version', version=__version__,
                        help="Print current stitchr version.")

    parser.add_argument('--cite', action=fxn.GetCitation, help="Print citation details.", nargs=0)

    return parser.parse_args()




def main():

    species = vars(args())['species'].upper()

    try:
        cmd = 'IMGTgeneDL -m stitchr -n -s ' + species
        subprocess.call(cmd, shell=True)
    except Exception:
        print("Failed to download sequences from IMGTgeneDL.")

    # Assuming that made the folder appropriately
    if species in os.listdir(os.getcwd()):
        # Delete older entries
        if species in os.listdir(fxn.data_dir):
            print("Note: a directory for this species exists: overwriting.")
            shutil.rmtree(os.path.join(fxn.data_dir, species))
        shutil.move(os.path.join(os.getcwd(), species), os.path.join(fxn.data_dir, species))

    else:
        raise IOError("Failed to produce the expected species\' output directory. ")
