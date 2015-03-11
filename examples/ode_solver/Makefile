# Make script to run and clean, etc.
# Please, be sure that you built the program.
#
# To run:
# $ make run
#
# To produce the plot AFTER the run:
# $ make plot
#
# To clean output files, plots, etc.
# $ make clean
#
# Vasily Eremenko <to@vsl.name>
# Linda Hlophe <...>
#
# 2015-03-11

SHELL = /bin/sh
MAKE = /usr/bin/make


###
SRCDIR = ../../src
PROG = ode
PLOTTER = vpot.py


###
.PHONY = all clean plot run


###
all :
	@echo Type \"make run\" to perform a run.


clean :
	rm -f k*.out
	rm -f vspot.out
	rm -f vkpot*.out
	rm -f vsep-srg.pdf


plot :
	python ${SRCDIR}/${PLOTTER}


run :
	time -p ${SRCDIR}/${PROG}
