
#CFLAGS=  -DCPLUSPLUS  -O2
CFLAGS= -g -G -DCPLUSPLUS  
INCL=  -I/data/disk03/huangk4/software/gsl-1.15/include
LIBS= /data/disk03/huangk4/software/gsl-1.15/lib/libgsl.a /data/disk03/huangk4/software/gsl-1.15/lib/libgslcblas.a  -lm 

all: wham_gpu.exe

wham.o: wham.cc wham.h
	nvcc $(INCL) $(CFLAGS) -c wham.cc 

wham_gpu.o: wham_gpu.cu wham.h
	nvcc $(INCL) $(CFLAGS) -c wham_gpu.cu

main.o: main.cc wham.h
	nvcc $(INCL) $(CFLAGS)  -c main.cc

wham_gpu.exe: wham.o main.o wham_gpu.o 
	nvcc   wham.o main.o wham_gpu.o  -o wham_gpu.exe $(LIBS)

clean: 
	rm *.exe *.o &>/dev/null
