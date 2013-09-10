
CFLAGS= -Wall -DCPLUSPLUS  -O2  
#CFLAGS= -Wall -DCPLUSPLUS  -g
INCL=  -I/home/huangk4/kun_gsl1.13/include
LIBS= /home/huangk4/kun_gsl1.13/lib/libgsl.a /home/huangk4/kun_gsl1.13/lib/libgslcblas.a  -lm 

wham.exe: wham.cc main.cc wham.h 
	g++ $(CFLAGS) $(INCL) wham.cc main.cc wham.h -o wham.exe $(LIBS)

clean: 
	rm *.exe *.o &>/dev/null
