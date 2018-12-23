#INCLUDES = syspara.h

bin_PROGRAMS = Voigt_atrial_2013
#Voigt_atrial_2013_CFLAGS = -O2 -Wall -axsse4.2 -mkl -parallel -par-report2 -Winline -inline-level=2 -ip -vec-report2 -opt-mem-layout-trans=3 -pg
Voigt_atrial_2013_CFLAGS = -O2 -Wall -axsse4.2 -mkl -parallel -par-report2 -vec-report2 -pg

Voigt_atrial_2013_INCLUDE = -I./include
Voigt_atrial_2013_LDADD = -L/opt/intel/mkl/lib/intel64 -liomp5 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core ./lib/xhplot.o

Voigt_atrial_2013_SOURCES = ./src/vGBd.c ./src/val_constants.c ./src/mem.c\
	./src/win_out.c ./src/eular.c ./src/runge.c ./src/input_para.c ./src/system.c ./src/currents.c\
	./src/eventloop.c ./src/exptable.c ./src/data_out.c ./src/syspara.h

