import numpy as np
import SimpleITK as sitk

L = 100 # mm total length in x, y, z
N = 100 # voxels per dimension
spacing = (L / N)

# coordinates along each axis
# CT/MHA store values at voxel centers so adding spacing/2

x = np.linspace(0, L, N, endpoint=False) + spacing/2
y = np.linspace(0, L, N, endpoint=False) + spacing/2
z = np.linspace(0, L, N, endpoint=False) + spacing/2

# Create 3D grids
# indexing ij keeps (z, y, x) in the correct CT-style order
zv, yv, xv = np.meshgrid(z, y, x, indexing="ij")

# CT (HU = 0 everywhere)
CT = np.zeros((N, N, N), dtype=np.int16)
img = sitk.GetImageFromArray(CT)
img.SetSpacing((spacing, spacing, spacing))
img.SetOrigin((x[0], y[0], z[0]))
sitk.WriteImage(img, "CT.mha")

# Bone slice
bone = np.zeros_like(CT)
bone[(xv >= 30) & (xv <= 70)] = 1

img2 = sitk.GetImageFromArray(bone)
img2.SetSpacing((spacing, spacing, spacing))
img2.SetOrigin((x[0], y[0], z[0]))
sitk.WriteImage(img2, "BoneSlice.mha")

print("Eureka! Eureka! It worked... probably!")
