###################################################################
# COPYRIGHT: EUMETSAT
#
# PRODUCED BY:
# Norwegian Meteorological Institute (DNMI)
# Research and Development Department
# P.O.BOX 43 - Blindern, N-0313 OSLO, NORWAY
#        
# This SW was developed by DNMI and DMI within the context of the
# Co-operation Agreement for the development of a pilot SAF on
# Ocean and Sea Ice.
###################################################################
#
# TYPE: 
# Makefile
# 
# PURPOSE:
# To create "accusnow", based on average_merge_SST
#
# AUTHOR:
# Steinar Eastwood, DNMI/FOU, 04.06.2001
# MODIFIED:
# Steinar Eastwood, met.no, 08.01.2003

# This include statement reads the module specification of required
# libraries and sets varaibles needed for LDFLAGS and CFLAGS. All required
# information for the module should be put in this file, even if not all
# parts of the module requires them. Specify the requirements by setting
# the adequate FLAGS below.

SHELL = /bin/sh

srcdir = .
prefix = /home/mariak/fmprojects/fmsnowcover/
exec_prefix = ${prefix}
bindir = ${exec_prefix}/bin
incdir = ${prefix}/include
libdir = ${exec_prefix}/lib

#CFLAGS = -g -O2 -DPACKAGE_NAME=\"METNO\ FMSNOWCOVER\" -DPACKAGE_TARNAME=\"metno-fmsnowcover\" -DPACKAGE_VERSION=\"0.1\" -DPACKAGE_STRING=\"METNO\ FMSNOWCOVER\ 0.1\" -DPACKAGE_BUGREPORT=\"o.godoy@met.no\" -DFMSNOWCOVER_HAVE_LIBPROJ=1 -DFMSNOWCOVER_HAVE_LIBHDF5=1 -DFMSNOWCOVER_HAVE_LIBMI=1 -DFMSNOWCOVER_HAVE_LIBUSENWP=1 -DFMSNOWCOVER_HAVE_LIBOSIHDF5=1 -DFMSNOWCOVER_HAVE_LIBTIFF=1 -DFMSNOWCOVER_HAVE_LIBFMUTIL=1 -DFMSNOWCOVER_HAVE_LIBFMIO=1 -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MALLOC=1  -I./

CFLAGS = -g -O2 -I./

CPPFLAGS =  -I/disk1/hdf5lib-1p6p//include -I/metno/local/lib -I/home/mariak/fmlocal/libs/libusenwp//include -I/home/mariak/fmlocal/libs/libosihdf5//include -I/home/mariak/fmlocal/libs/libfmutil//include -I/home/mariak/fmlocal/libs/libfmio//include

LDFLAGS =  -L/disk1/hdf5lib-1p6p//lib -L/metno/local/lib -L/home/mariak/fmlocal/libs/libusenwp//lib -L/home/mariak/fmlocal/libs/libosihdf5//lib -L/home/mariak/fmlocal/libs/libfmutil//lib -L/home/mariak/fmlocal/libs/libfmio//lib

LIBS = -Wl,-rpath=/home/mariak/fmlocal/libs/libfmio//lib -lfmio -Wl,-rpath=/home/mariak/fmlocal/libs/libfmutil//lib -lfmutil -ltiff -Wl,-rpath=/home/mariak/fmlocal/libs/libosihdf5//lib -losihdf5 -Wl,-rpath=/home/mariak/fmlocal/libs/libusenwp//lib -lusenwp -Wl,-rpath=/metno/local/lib -lmi -Wl,-rpath=/disk1/hdf5lib-1p6p//lib -lhdf5 -lproj  -lg2c


HEADER_FILES = \
  accusnow.h 

SRC_FILES = \
  accusnow.c \
  store_snow.c \
  func_accusnow.c 
  
BINFILE = accusnow

OBJ_FILES := $(SRC_FILES:.c=.o)

$(BINFILE): $(OBJ_FILES) 
	$(CC) $(CFLAGS) -o $(BINFILE) $^ $(LDFLAGS) $(LIBS)

$(OBJ_FILES): $(HEADER_FILES)



# Cleaning and installing commands
clean:
	-rm -f $(OBJS)

distclean:
	-rm -f $(OBJS)
	-rm -f $(RUNFILE)
	-rm -f ../bin/$(RUNFILE)

install:
	install -d ../bin
ifdef RUNFILE
	install $(RUNFILE) ../bin
endif
ifdef JOBFILES
	install $(JOBFILES) ../job
endif
