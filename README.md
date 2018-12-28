# kirchhoff_verlet

This is a small code that simulates the packing of long slender rods into spherical cavities, as used, e.g. in Stoop et al, Phys. Rev. Lett. 106, 214102 (2011). The slender rod is discretized into a set of N segments and moments and forces of each segment are calculated similarly to J. Spillmann et al.,  Proc. ACM SIGGRAPH/Eurographics (2007).

This code requires the MILI header libraries, available at https://bitbucket.org/fudepan/mili/wiki/Home You will need to adapt the Makefile to include/link MILI appropriately.

