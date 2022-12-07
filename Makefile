include make.inc

TARGET = SSBE 

OBJS = \
src/math/math_constants.o \
src/math/phys_constants.o \
src/parallel/communication.o \
src/io/salmon_file.o \
src/common/structures.o \
src/common/pack_unpack.o \
src/util.o \
src/sbe_gs.o \
src/sbe_bloch_solver.o \
src/test.o \
src/input_parameter.o \
src/rt/em_field.o \
src/maxwell/fdtd_weyl_gauge.o \
src/realtime.o \
src/multiscale.o \
src/main.o

$(TARGET): $(OBJS)
	$(FC) -o $@ $^ $(FLAGS) $(LIBS)

.SUFFIXES: .f90

%.o: %.f90
	$(FC) -c -o $@ $^ $(FLAGS)

.PHONY: all clean 

all: $(TARGET)

clean:
	rm $(TARGET) $(OBJS) *.mod

#input_parameter.f90: input_parameter.py
#	python input_parameter.py > $@
