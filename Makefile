all: main

main:
	g++ cleaning_corsika_auto.cpp -o cleaning -Wall -Wfloat-equal -Warray-bounds -Wdiv-by-zero -Wextra -Wconversion
	g++ find_corsica_stereo_auto.cxx -o find_stereo -Wall -Wfloat-equal -Warray-bounds -Wdiv-by-zero -Wextra -Wconversion
	g++ trigger_iact234.cpp -o trigger_iact234 -lMathMore -fopenmp -I `root-config --incdir` `root-config --cflags` `root-config --libs`
	g++ trigger_iact1.cpp -o trigger_iact1 -lMathMore -fopenmp -I `root-config --incdir` `root-config --cflags` `root-config --libs`
clean:
	rm trigger_iact234 && trigger_iact1 && find_stereo && cleaning
