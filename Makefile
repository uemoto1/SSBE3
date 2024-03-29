include make.inc

TARGET = SSBE 

OBJS = \
src/math/math_constants.o \
src/math/phys_constants.o \
src/parallel/communication.o \
src/io/salmon_file.o \
src/input_parameter.o \
src/common/structures.o \
src/common/pack_unpack.o \
src/util.o \
src/sbe_gs.o \
src/sbe_bloch_solver.o \
src/test.o \
src/rt/em_field.o \
src/maxwell/fdtd_weyl_gauge.o \
src/realtime.o \
src/multiscale.o \
src/main.o

$(TARGET): $(OBJS)
	$(FC) -o $@ $^ $(FLAGS) $(LIBS)

%.o: %.f90
	$(FC) -c -o $@ $^ $(FLAGS)

.PHONY: all clean 

all: $(TARGET)

clean:
	rm $(TARGET) $(OBJS) src/*.mod
