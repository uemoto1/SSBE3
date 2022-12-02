include make.inc

TARGET = SSBE 

OBJS = \
math/math_constants.o \
math/phys_constants.o \
parallel/communication.o \
io/salmon_file.o \
common/structures.o \
common/pack_unpack.o \
util.o \
sbe_gs.o \
sbe_bloch_solver.o \
test.o \
input_parameter.o \
rt/em_field.o \
maxwell/fdtd_weyl_gauge.o \
realtime.o \
multiscale.o \
main.o

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
