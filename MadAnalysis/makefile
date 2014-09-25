ifeq ($(origin FFLAGS),undefined)
  FFLAGS= -O -w
  #FFLAGS+= -g -fbounds-check -ffpe-trap=invalid,zero,overflow,underflow,denormal -Wall
endif
FFLAGS+= -ffixed-line-length-132

# Set FC unless it's defined by an environment variable
ifeq ($(origin FC),default)
  FC=g77 
endif

PLOT    =  plot_events.o dbook.o rw_events.o kin_func.o set_cuts.o

.f.o: ; $(F77) $(FFLAGS) -c $*.f

plot_events: $(PLOT)
	$(F77) $(FFLAGS) -o plot_events $(PLOT)

kin_func.o: info.inc kin_func.f
plot_events.o: info.inc plot_events.f
set_cuts.o: info.inc set_cuts.f
dbook.o: dbook.inc dbook.f

clean:
	rm -f *.o
	rm -f plot_events
