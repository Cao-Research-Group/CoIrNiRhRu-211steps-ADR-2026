(i) runMonteCanonical.java: The code to run canonical MMC simulations;

(ii) DecorationMSymCanonicalMidFreezeManager.java: MMC manager, with middle 1 layer fixed and other layers canonical, for mirror symmetric slab;

(iii) DecorationMSymCanonicalManager.java: MMC manager, canonical, for mirror symmetric slab;

(iv) CombinedRecorder.java: Combined recorder integrating site composition, *N adsorption energy, and TOF prediction.

(v) Sym_equimolar-prim=80.vasp: An initial structure for example, which is an equimolar 15 x 16 slab.

(vi) PRIM.vasp: the PRIM file for CoIrNiRhRu fcc-(211) CE.

In our settings, The annealing schedule was set to cool from 2000 K to 1000 K in ¦¤T = 20 K increments. 
At each temperature, the equilibration length of the MMC simulations was set to 4500 times of the total sites;
the recording length was set to 500 times of the total sites. 
For our 15¡Á16 slab (1840 sites), this corresponds to 8,280,000 equilibration steps and 920,000 recording steps per temperature.
