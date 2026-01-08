import sys
sys.path.append('..')
import os

from opentps.core.data.plan import ProtonPlan, RTPlan
from opentps.core.io.scannerReader import readScanner
from opentps.core.io import mcsquareIO
from opentps.core.processing.doseCalculation.doseCalculationConfig import DoseCalculationConfig
from opentps.core.processing.doseCalculation.protons.mcsquareDoseCalculator import MCsquareDoseCalculator
from opentps.core.data.plan._planProtonBeam import PlanProtonBeam
from opentps.core.data.plan._planProtonLayer import PlanProtonLayer
from opentps.core.data.images._ctImage import CTImage
import numpy as np
import matplotlib.pyplot as plt

ctCalibration = readScanner(DoseCalculationConfig().scannerFolder)
bdl = mcsquareIO.readBDL(DoseCalculationConfig().bdlFile)

# Loading a copy of bdl file to prevent any modification of 
# ctCalibration = readScanner(DoseCalculationConfig().scannerFolder)
# bdl_path = "/home/amit/opentps/opentps_core/opentps/core/processing/doseCalculation/protons/MCsquare/BDL/BDL_fred_like_pbwater.txt"
# bdl = mcsquareIO.readBDL(bdl_path)
ctSize = 300
#bdl.nozzle_isocenter = 400 # making at same as that of FRED to match the beam model

ct = CTImage()
ct.name = 'CT'
huAir = -1024.
huWater = ctCalibration.convertRSP2HU(1.)
data = huWater * np.ones((ctSize, ctSize, ctSize))
#data[:, :50, :] = huAir
ct.imageArray = data

ct.spacing = (1.0, 1.0, 1.0)
voxel_spacing = float(ct.spacing[1])

plan = ProtonPlan()
plan.appendBeam(PlanProtonBeam())
plan.beams[0].gantryAngle = 0
plan.beams[0].appendLayer(PlanProtonLayer(150)) # 150 MeV
plan[0].layers[0].appendSpot([0], [0], [1])

# Why do you do this ctSize*voxel_spacing//2.0 and not just ctSize//2.0?
# total size = N x spacing
plan.beams[0].isocenterPosition = [ctSize*voxel_spacing//2.0, ctSize*voxel_spacing//2.0, ctSize*voxel_spacing//2.0]


#Making the beam as pencil beam with spot size = 0
# to match the beam parameters (twiss) that of fred
#sigma_x = 1.69 # mm
sigma_spatial = 0.02
bdl.spotSize1x = [sigma_spatial] * len(bdl.nominalEnergy)
bdl.spotSize1y = [sigma_spatial] * len(bdl.nominalEnergy)
bdl.spotSize2x = [sigma_spatial] * len(bdl.nominalEnergy)
bdl.spotSize2y = [sigma_spatial] * len(bdl.nominalEnergy)

# to match the twiss parameter, because this can't be set to zero
#sigma_div = 0.0125 # rad, dimensionless
sigma_div = 0.2e-3
bdl.divergence1x = [sigma_div] * len(bdl.nominalEnergy)
bdl.divergence1y = [sigma_div] * len(bdl.nominalEnergy)
bdl.divergence2x = [sigma_div] * len(bdl.nominalEnergy)
bdl.divergence2y = [sigma_div] * len(bdl.nominalEnergy)

# For making EFWHM = 0 MeV
bdl.energySpread = []
for e in bdl.nominalEnergy:
    spread = (2.12 * 100)/e #Mono energetic beam
    bdl.energySpread.append(spread)


mc2 = MCsquareDoseCalculator()
mc2.beamModel = bdl
mc2.ctCalibration = ctCalibration
mc2.nbPrimaries = 1e6
dose = mc2.computeDose(ct, plan)

# print(plan.beams[0].isocenterPosition)
# print(plan.beams[0].mcsquareIsocenter)

# print(dir(bdl))
# print(bdl.nozzle_isocenter)

#plt.imshow(dose.imageArray[75, :, :])
#plt.show()
#print(dir(bdl))

# output_dir = "/home/amit/Plots"
# os.makedirs(output_dir, exist_ok=True)

# output_path = os.path.join(output_dir, "dose_slice1.png")
# plt.imshow(dose.imageArray[150, :, :])
# plt.title("Dose slice at _ [150, :, :]")
# plt.colorbar(label="Dose (Gy)")
# plt.savefig(output_path, dpi=300, bbox_inches="tight")
# plt.close()
# print(f"Saved plot at: {output_path}")

# output_path = os.path.join(output_dir, "dose_slice2.png")
# plt.imshow(dose.imageArray[:, 150, :])
# plt.title("Dose slice at _ :, 150, :]")
# plt.colorbar(label="Dose (Gy)")
# plt.savefig(output_path, dpi=300, bbox_inches="tight")
# plt.close()
# print(f"Saved plot at: {output_path}")

# output_path = os.path.join(output_dir, "dose_slice3.png")
# plt.imshow(dose.imageArray[:, :, 150])
# plt.title("Dose slice at _ [:, :, 150]")
# plt.colorbar(label="Dose (Gy) ")
# plt.savefig(output_path, dpi=300, bbox_inches="tight")
# plt.close()
# print(f"Saved plot at: {output_path}")