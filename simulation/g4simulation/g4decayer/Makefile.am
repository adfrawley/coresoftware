AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libg4decayer.la

INCLUDES = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I${G4_MAIN}/include \
  -I${G4_MAIN}/include/Geant4

libg4decayer_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OPT_SPHENIX)/lib \
  `geant4-config --libs`

libg4decayer_la_LIBADD = \
  -lPythia6

libg4decayer_la_SOURCES = \
  pythia6_common_address.c \
  Pythia6.cc \
  G4Pythia6Decayer.cc \
  G4Pythia6DecayerMessenger.cc \
  P6DExtDecayerPhysics.cc \
  P6DPhysicsList.cc

pkginclude_HEADERS = \
  EDecayType.hh \
  G4Pythia6Decayer.hh \
  G4Pythia6DecayerMessenger.hh \
  P6DExtDecayerPhysics.hh \
  P6DPhysicsList.hh

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.C

testexternals_LDADD = \
  libg4decayer.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict*
