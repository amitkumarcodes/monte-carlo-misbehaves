import sys
sys.path.append('..')
import os
import SimpleITK as sitk

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
# ctSize = 100


# ct = CTImage()
# ct.name = 'CT'
# huAir = -1024.
# huFakeBone = 1300
# huWater = ctCalibration.convertRSP2HU(1.)
# data = huWater * np.ones((ctSize, ctSize, ctSize))

# ct.origin = (0, 0, 0)
# ct.spacing = (1, 1, 1)
# voxel_spacing = float(ct.spacing[1])

# # Designing the Heterogenous Phantom
# L = ctSize * voxel_spacing
# x = np.linspace(-L/2, L/2, ctSize, endpoint = False) + voxel_spacing//2
# y = np.linspace(-L/2, L/2, ctSize, endpoint = False) + voxel_spacing//2
# z = np.linspace(0, L, ctSize, endpoint = False) + voxel_spacing//2

# # create 3D grids
# # indexing ij keeps (z, y, x) in the correct CT-style order
# zv, yv, xv = np.meshgrid(z, y, x, indexing="ij")

# # 3 cm thick bone slab 
# bone_mask = (zv > 30) & (zv < 70)

# data[bone_mask] = huFakeBone
# ct.imageArray = data

############# CT images from FRED ##############
ct_sitk = sitk.ReadImage("CT.mha")
ct_arr = sitk.GetArrayFromImage(ct_sitk)
bone_sitk = sitk.ReadImage("BoneSlice.mha")
bone_arr = sitk.GetArrayFromImage(bone_sitk)

ct = CTImage()
ct.name = "CT"
ct.imageArray = ct_arr
ct.spacing = ct_sitk.GetSpacing()
ct.origin = ct_sitk.GetOrigin()

Bone = 1300.0
bone_mask = bone_arr > 0.5
ct.imageArray[bone_mask] = Bone

nz, ny, nx = ct.imageArray.shape
dz, dy, dx = ct.spacing
ox, oy, oz = ct.origin

iso_x = ox + dx * nx / 2.0
iso_y = oy + dy * ny / 2.0
iso_z = oz + dz * nz / 2.0

plan = ProtonPlan()
plan.appendBeam(PlanProtonBeam())
plan.beams[0].gantryAngle = 0
plan.beams[0].appendLayer(PlanProtonLayer(100)) # 100 MeV

#plan[0].layers[0].appendSpot([0], [0], [1])

#### Proton confetti
layer = plan.beams[0].layers[0]

coords = np.linspace(-40, 40, 40)
X, Y = np.meshgrid(coords, coords)

x = X.ravel()
y = Y.ravel()
mu = np.ones_like(x)

for xi, yi, mui in zip(x, y, mu):
    plan.beams[0].layers[0].appendSpot([xi], [yi], [mui]) ###changed this one


# Why do you do this ctSize*voxel_spacing//2.0 and not just ctSize//2.0?
# total size = N x spacing
#plan.beams[0].isocenterPosition = [ctSize*voxel_spacing//2.0, ctSize*voxel_spacing//2.0, ctSize*voxel_spacing//2.0]
plan.beams[0].isocenterPosition = [iso_x, iso_y, iso_z]

#Making the beam as pencil beam with spot size = 0
tiny = 2.0
bdl.spotSize1x = [tiny] * len(bdl.nominalEnergy)
bdl.spotSize1y = [tiny] * len(bdl.nominalEnergy)
bdl.spotSize2x = [tiny] * len(bdl.nominalEnergy)
bdl.spotSize2y = [tiny] * len(bdl.nominalEnergy)

tiny_div = 1e-9
bdl.divergence1x = [tiny_div] * len(bdl.nominalEnergy)
bdl.divergence1y = [tiny_div] * len(bdl.nominalEnergy)
bdl.divergence2x = [tiny_div] * len(bdl.nominalEnergy)
bdl.divergence2y = [tiny_div] * len(bdl.nominalEnergy)

# For making EFWHM = 5 MeV
bdl.energySpread = []
for e in bdl.nominalEnergy:
    spread = (2.12 * 100)/e
    bdl.energySpread.append(spread)


mc2 = MCsquareDoseCalculator()
mc2.beamModel = bdl
mc2.ctCalibration = ctCalibration
mc2.nbPrimaries = 5e7
dose = mc2.computeDose(ct, plan)

print(mc2.nbPrimaries)
# print(plan.beams[0].isocenterPosition)
# print(plan.beams[0].mcsquareIsocenter)

#print(dir(bdl))

# plt.imshow(dose.imageArray[75, :, :])
# plt.show()
# print(dir(bdl))

# output_dir = "/home/amit/Plots"
# os.makedirs(output_dir, exist_ok=True)

# output_path = os.path.join(output_dir, "dose_slice1.png")
# plt.imshow(dose.imageArray[50, :, :])
# plt.title("Dose slice at _ [50, :, :]")
# plt.colorbar(label="Dose (Gy)")
# plt.savefig(output_path, dpi=300, bbox_inches="tight")
# plt.close()
# print(f"Saved plot at: {output_path}")

# output_path = os.path.join(output_dir, "dose_slice2.png")
# plt.imshow(dose.imageArray[:, 50, :])
# plt.title("Dose slice at _[:, 50, :]")
# plt.colorbar(label="Dose (Gy)")
# plt.savefig(output_path, dpi=300, bbox_inches="tight")
# plt.close()
# print(f"Saved plot at: {output_path}")

# output_path = os.path.join(output_dir, "dose_slice3.png")
# plt.imshow(dose.imageArray[:, :, 50])
# plt.title("Dose slice at _ [:, :, 50]")
# plt.colorbar(label="Dose (Gy)")
# plt.savefig(output_path, dpi=300, bbox_inches="tight")
# plt.close()
# print(f"Saved plot at: {output_path}")

# output_path = os.path.join(output_dir, "idd.png")
# idd = dose_gy.imageArray.sum(axis=(0,2))
# plt.plot(idd)
# plt.savefig(output_path, dpi=300, bbox_inches="tight")
# plt.close()
# print(f"Saved plot at: {output_path}")