FFOPTS = -g -ftrap=common -c
#FFOPTS =  -c -C
#FFOPTS = -c -g -C -f77=%all,logical  -onetrip -Bstatic 

.SUFFIXES: .o .f .f90

OBJS   = newtonr.o LippSchw.o potential.o
OBJS  += gauleg.o 

#Change the location to the location of lapack and blas libraries
LIBDIR = /home/sudhanva/Documents/lib

#----------------------------------------------------

all : lipp

lipp : $(OBJS)
	f95 -o lipp $(OBJS) -L$(LIBDIR)  -llapack_LINUX -lblas_LINUX

%.o : %.f90
	f95 $(FFOPTS) $^

clean:
	@rm -rf *o *~ lipp
	@echo cleaned

