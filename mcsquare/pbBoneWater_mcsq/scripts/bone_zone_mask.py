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

ctSize = 100

ct = CTImage()
ct.name = 'CT'
huAir = -1024
huBone = 1300
huWater = ctCalibration.convertRSP2HU(1.)
data = huWater * np.ones((ctSize, ctSize, ctSize))
data[:, 30:70, :] = huBone
ct.imageArray = data
ct.origin = (0, 0, 0)
ct.spacing = (1.0, 1.0, 1.0)
voxel_spacing = float(ct.spacing[1])

plan = ProtonPlan()
plan.appendBeam(PlanProtonBeam())
plan.beams[0].gantryAngle = 0
plan.beams[0].appendLayer(PlanProtonLayer(100))
plan[0].layers[0].appendSpot([0], [0], [1])

plan.beams[0].isocenterPosition = [ctSize*voxel_spacing//2.0, ctSize*voxel_spacing//2.0, ctSize*voxel_spacing//2.0]

tiny = 1e-3
bdl.spotSize1x = [tiny] * len(bdl.nominalEnergy)
bdl.spotSize1y = [tiny] * len(bdl.nominalEnergy)
bdl.spotSize1z = [tiny] * len(bdl.nominalEnergy)

bdl.divergence1x = [tiny] * len(bdl.nominalEnergy)
bdl.divergence1y = [tiny] * len(bdl.nominalEnergy)
bdl.divergence2x = [tiny] * len(bdl.nominalEnergy)
bdl.divergence2y = [tiny] * len(bdl.nominalEnergy)

bdl.energySpread = []
for e in bdl.nominalEnergy:
    spread = (2.12 * 100)/e
    bdl.energySpread.append(spread)

mc2 = MCsquareDoseCalculator()
mc2.beamModel = bdl
mc2.ctCalibration = ctCalibration
mc2.nbPrimaries = 1e6
dose = mc2.computeDose(ct, plan)

print(plan.beams[0].isocenterPosition)
print(plan.beams[0].mcsquareIsocenter)

print(huWater)
print(huBone)