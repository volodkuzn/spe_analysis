import os

import matplotlib.pyplot as plt

import spe
import spe.spectrum


os.chdir("../")
sp = spe.spectrum.Spectrum2D("./Sasha/8-240_polarization/11lum_-8_0_0d025_8002_.SPE",
                             spe.spectrum.Axis("Field", "T", 0, -8), max_intensity=25000)
sp.remove_cosmic_rays()
sp.plot_all(figsize=(20, 10))
plt.show()