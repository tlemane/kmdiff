override CFLAGS += -I../include -I/usr/include/openblas
override LDLIBS += -lopenblas -lgsl -lopenblas -lm -lpthread 
override LDLIBS += -llapacke

override CFLAGS += -O2

ifeq ($(DEBUG), 1)
	CFLAGS += -g # enable debugging
endif

ifeq ($(PROFILING), 1)
	CFLAGS += -pg # enable profiling
endif

ND=nicksrc
ED=eigensrc
KD=ksrc

NLIB = $(ND)/libnick.a

# ----- phony targets
.PHONY: all clean clobber install

EXE = baseprog convertf mergeit pca \
	$(ED)/pcatoy $(ED)/smartrel $(ED)/smarteigenstrat \
	$(ED)/twstats $(ED)/eigenstrat $(ED)/eigenstratQTL $(ED)/smartpca \
	smartpca 

all:    $(EXE)

install:  all
	mv $(EXE) ../bin

clobber:
	rm -f *.o */*.o */*.a $(ND)/*.o $(ED)/*.o $(KD)/*.o
	rm -f $(EXE)
	cd ../bin/ ; rm -f $(notdir $(EXE)) ; cd ../src

clean:
	rm -f *.o core core.* *.o $(ND)/*.o $(ED)/*.o $(KD)/*.o $(EXE)

# ----- build nicksrc/libnick.a
$(NLIB):
	$(MAKE) -C $(ND)

baseprog:	baseprog.o mcio.o egsubs.o admutils.o h2d.o $(ED)/exclude.o $(NLIB)

convertf:	convertf.o mcio.o egsubs.o admutils.o h2d.o $(ED)/exclude.o $(NLIB)

mergeit:	mergeit.o mcio.o admutils.o $(NLIB)

pca:	pca.o $(ED)/eigsubs.o eigensrc/eigx.o $(NLIB)

$(ED)/pcatoy: 	$(ED)/pcatoy.o eigensrc/eigsubs.o eigensrc/eigx.o $(NLIB)

$(ED)/smartrel:	$(ED)/smartrel.o twsubs.o mcio.o qpsubs.o admutils.o egsubs.o regsubs.o \
		$(ED)/eigsubs.o $(ED)/eigx.o $(ED)/smartsubs.o $(NLIB)    

$(ED)/smarteigenstrat:	$(ED)/smarteigenstrat.o mcio.o admutils.o $(NLIB)

$(ED)/twstats:	$(ED)/twstats.o $(NLIB) 

$(ED)/smartpca:	$(ED)/smartpca.o $(ED)/eigsubs.o $(ED)/exclude.o $(ED)/smartsubs.o $(ED)/eigx.o \
		mcio.o qpsubs.o admutils.o egsubs.o regsubs.o gval.o \
		$(NLIB) \
          	$(KD)/kjg_fpca.o $(KD)/kjg_gsl.o twsubs.o

smartpca:   $(ED)/smartpca
	cp $(ED)/smartpca smartpca
