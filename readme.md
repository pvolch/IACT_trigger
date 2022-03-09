This repository contains several scripts.

1) readbin scripts are designed to introduce a telescope trigger and night sky background into simulated events recorded by telescopes. The data acquisition system of the first telescope is different from the others, so two scripts are present.

To run, you need to specify the path and name of the source file and the path to the out folder in the script itself. A description of the input binary file can be found on the GridMS server. The specification of the output files is also located there.

2) cleaning_corsika_auto.cpp completes the readbin output files with cleaning with the thresholds specified in the "clean.param" file. This file contains the telescopes to be cleaned and other specifications.

3)find_corsica_stereo_auto.cxx - generates a common installation joint events file. Uses the specifications specified in "clean.param".
