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

ctSize = 300
bdl.nozzle_isocenter = 400 # making at same as that of FRED to match the beam model

ct = CTImage()
ct.name = 'CT'
huAir = -1024
data = huAir * np.ones((ctSize, ctSize, ctSize))
ct.imageArray = data

ct.spacing = (1.0, 1.0, 1.0)
voxel_spacing = float(ct.spacing[1])

plan = ProtonPlan()
plan.appendBeam(PlanProtonBeam())
plan.beams[0].gantryAngle = 0
plan.beams[0].appendLayer(PlanProtonLayer(1))
plan[0].layers[0].appendSpot([0], [0], [1])

#total size = N * spacing
plan.beams[0].isocenterPosition = [ctSize*voxel_spacing//2.0, ctSize*voxel_spacing//2.0, ctSize*voxel_spacing//2.0]

#making the beam as spot size = 0 almost
sigma_spatial = 1e-9
bdl.spotSize1x = [sigma_spatial] * len(bdl.nominalEnergy)
bdl.spotSize1y = [sigma_spatial] * len(bdl.nominalEnergy)
bdl.spotSize2x = [sigma_spatial] * len(bdl.nominalEnergy)
bdl.spotSize2y = [sigma_spatial] * len(bdl.nominalEnergy)

sigma_div = 1e-9
bdl.divergence1x = [sigma_div] * len(bdl.nominalEnergy)
bdl.divergence1y = [sigma_div] * len(bdl.nominalEnergy)
bdl.divergence2x = [sigma_div] * len(bdl.nominalEnergy)
bdl.divergence2y = [sigma_div] * len(bdl.nominalEnergy)

bdl.energySpread = []
for e in bdl.nominalEnergy:
    spread = 0 * (2.12 * 100)/e
    bdl.energySpread.append(spread)

mc2 = MCsquareDoseCalculator()
mc2.beamModel = bdl
mc2.ctCalibration = ctCalibration
mc2.nbPrimaries = 1e6
dose = mc2.computeDose(ct, plan)