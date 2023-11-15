F90 = mpiifort
Exe = bicycle

# Files

Files = Bicycle_2D.f

Obs =  $(Files:.f=.o)

Libs = -lfftw3

Incs = 

Opt = -O3 -mcmodel=medium -shared-intel -extend-source 132 

$(Exe): $(Obs)
	$(F90) $(Opt)$ -o $(Exe) $(Obs) $(Libs)

%.o: %.f
	$(F90) -c $(Opt) $(Incs) $<

clean:
	rm -f *.o *.out
