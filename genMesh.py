import sys, os
import numpy as np
import argparse
from matplotlib import pyplot as plt
from pyhyp import pyHyp
from cgnsutilities import cgnsutilities
from prefoil import Airfoil, sampling

# Add the path so we can import constants directly
sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
from SETUP import constants

parser = argparse.ArgumentParser()
parser.add_argument(
    "--gridFamily",
    help="Family of grids",
    default="A",
    choices=["A", "B"],
)
parser.add_argument(
    "--input",
    help="Input airfoil file",
    type=str,
    default="OAT15A.dat",
)
parser.add_argument(
    "--geometry",
    help="Input geometry configuration.",
    default="G2",
    choices=["G1", "G2"],
    type=str,
)
parser.add_argument("--debug", help="Check point sampling", action="store_true")
args = parser.parse_args()

geometry = {
    "G1": {"chord": 0.23},
    "G2": {"chord": 1.0},#"G2": {"chord": 1.0},
}

chord = geometry[args.geometry]["chord"]
if args.gridFamily == "A":
    nPts = 145#293
    N = 85
    nTEPts = 11
    s0 = 1e-6#5e-7
elif args.gridFamily == "B":
    nPts = 897
    N = 97
    nTEPts = 15
    s0 = 1e-6

if not args.input:
    raise Warning("Provide input coordinates")

# Strip first col from the data
coords = np.loadtxt(args.input)
coords = coords[:, 1:]

airfoilName = args.input.replace(".dat", "")
airfoilName += f"_{args.geometry}_{args.gridFamily}"

# Sample points using prefoil
airfoil = Airfoil(coords, normalize=True)
airfoil.scale(factor=chord)
data = airfoil.getSampledPts(nPts, spacingFunc=sampling.conical, func_args={"coeff": 1}, nTEPts=nTEPts)

if args.debug:  # checking out the airfoil and point distribution
    fig1 = airfoil.plot()
    fig1.suptitle(f"{airfoilName} with selected coordinates sampling")
    plt.show()

# Write the coords
airfoil.writeCoords(airfoilName)
plot3DFileName = f"{airfoilName}.xyz"

options = {
    # ---------------------------
    #        Input Parameters
    # ---------------------------
    "inputFile": "plot3DFileName", #"RAE2822_spanwise21_scaled.xyz",
    "unattachedEdgesAreSymmetry": True,
    "outerFaceBC": "farfield",
    "autoConnect": True,
    "BC": {1: {"jLow": "zSymm", "jHigh": "zSymm"}},
    "families": "wall",
    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    "N": N,
    "s0": s0*chord,
    "marchDist": 100*chord,
    #"nConstantStart": 10,
    
}

cgnsFileName = f"{airfoilName}_L0.cgns"
hyp = pyHyp(options=options)
hyp.run()
hyp.writeCGNS(cgnsFileName)

'''
# Coarsen the mesh
grid = cgnsutilities.readGrid(cgnsFileName)
grid.coarsen()
grid.writeToCGNS(f"{airfoilName}_L1.cgns")
grid.coarsen()
grid.writeToCGNS(f"{airfoilName}_L2.cgns")
'''
