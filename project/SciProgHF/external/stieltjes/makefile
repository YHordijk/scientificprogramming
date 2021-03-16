CC=gcc
F90=gfortran
F77=gfortran
#
Dd0=-DDEBUG # general debugging
Dd1=-DDBGCLASS # cross-check the configs classification
Dd2=-DFULLHAM # calculate the full spectrum using full matrix
Dd3=-DDBGSTIELTJES # still more output from Stieltjes
Dd4=-DDBGPARTIAL # debug info about the partial widths calculation
#
D1=-DEXTOUT
D2=-DGERES # use e_init for E in <in|(H-E)|fin_j> instead of e_fin(j)
D3=-DESHIFT # shift e_fin by e_init and use offset
D4=-DFULLINEXP # store weights of leading configurations for complete IN spectrum
#
DEFS=$(D2) $(D3) -DADAPTDSYEV
TARG=stieltjes.x

OBJECTS= herm_spline.o stieltjes_ql_diag.o stieltjesmod.o input_routines.o driver.o

$(TARG): $(OBJECTS) makefile
	$(F90) $(OBJECTS) -o $(TARG)

herm_spline.o: herm_spline.c makefile
stieltjes_ql_diag.o: stieltjes_ql_diag.f90 makefile
stieltjesmod.o: herm_spline.o stieltjes_ql_diag.o stieltjesmod.f90 makefile
input_routines.o: input_routines.f90 makefile
driver.o: input_routines.o stieltjesmod.o driver.f90 makefile

%.o: %.f90
	$(F90) $(OPTS) -c $<

%.o: %.f
	$(F77) $(OPTS) -c $<

%.o: %.c
	$(CC) $(COPTS)  -c $<

clean:
	rm -f *.o *.mod *__genmod.f90

fullclean:
	rm -f *.o *.mod *__genmod.f90 $(TARG)
