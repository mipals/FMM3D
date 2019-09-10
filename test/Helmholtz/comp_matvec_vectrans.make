PROJECT = int2

HOST = macosx
HOST = linux-gfortran
HOST = linux-ifort
HOST = linux-gfortran-openmp

ifeq ($(HOST),macosx)
FC = gfortran -c -w
FFLAGS = -O3
FLINK = gfortran -w -o $(PROJECT) -framework accelerate
endif

ifeq ($(HOST),linux-gfortran)
FC = gfortran -c 
FFLAGS = -O3 -march=native -funroll-loops -fPIC -ftree-vectorizer-verbose=2
FLINK = gfortran -o $(PROJECT) -lopenblas
endif

ifeq ($(HOST),linux-ifort)
FC = ifort -c 
FFLAGS = -xHost -O3 -g -qopenmp 
FLINK = ifort -o $(PROJECT) -mkl -qopenmp
endif

ifeq ($(HOST),linux-gfortran-openmp)
FC = gfortran 
FFLAGS = -O3 -march=native -fPIC -funroll-loops  -c --openmp
FLINK = gfortran -o $(PROJECT) --openmp -lopenblas
endif



COM = ../../src/Common
HELM = ../../src/Helmholtz


.PHONY: all clean list


SOURCES =  comp_matvec_vectrans.f \
  $(COM)/hkrand.f \
  $(COM)/dlaran.f \
  $(COM)/prini.f \
  $(COM)/rotviarecur.f \
  $(HELM)/h3dtrans.f \
  $(COM)/dfft.f \
  $(COM)/rotproj.f \
  $(HELM)/projections.f \
  $(HELM)/h3dcommon.f \
  $(COM)/yrecursion.f \
  $(COM)/fmmcommon.f \
  $(COM)/besseljs3d.f \
  $(COM)/legeexps.f \
  $(HELM)/h3dterms.f \
  $(HELM)/helmrouts3d.f 


OBJECTS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCES)))
#

%.o : %.f
	$(FC) $(FFLAGS) $< -o $@

%.o : %.f90
	$(FC) $(FFLAGS) $< -o $@

all: $(OBJECTS)
	rm -f $(PROJECT)
	$(FLINK) $(OBJECTS) $(FEND)
	./$(PROJECT)

clean:
	rm -f $(OBJECTS)
	rm -f $(PROJECT)

list: $(SOURCES)
	$(warning Requires:  $^)





  
