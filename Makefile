########################################################################################################
# Intel fortran Compilateur
# OPTION A EVITER AVEC IFORT: -fast
#F90      = ifort
# Options de compilation [mode normal] 
#F90FLAGS = -zero -autodouble -i8 -r8 -O2 -fPIC -ipo -inline all -implicitnone -fp-model precise
#F90FLAGS = -warn unused -check all -zero -autodouble -i8 -r8 -O2 -fPIC -ipo -inline all -openmp -implicitnone -fp-model precise
#F90FLAGS = -check all -zero -autodouble -i8 -r8 -O2 -fPIC -ipo -inline all -implicitnone -fp-model precise
# Options de compilation [mode paranoiaque]
#F90FLAGS = -r8 -i8 -autodouble -O0 -fno-alias -fno-fnalias -check all -implicitnone -warn all  -fp-stack-check  \
           -fPIC -zero -ftrapuv -heap-arrays -gen-interface 
#LDFLAGS  = -fPIC
#LDFLAGS  = -openmp -fPIC
########################################################################################################
# GNU fortran compileur
F90      = gfortran
# Options de compilation [mode normal] 
F90FLAGS = -O2 -fdefault-integer-8 -fdefault-real-8 -fimplicit-none -fPIC -ftree-vectorize 
#F90FLAGS = -O2 -fdefault-integer-8 -fdefault-real-8 -fPIC
#F90FLAGS = -O2 -ftree-vectorize -ftree-vectorizer-verbose=2 -fbounds-check -fdefault-integer-8  -fimplicit-none -Wconversion -fPIC \
	   -Wsurprising -Wno-tabs -Wunderflow -Wintrinsic-shadow -Werror -Wunused-parameter
# Options de compilation [mode paranoiaque]
#F90FLAGS = -finit-local-zero -fbounds-check -fimplicit-none \
           -pedantic -pedantic-errors -Wall -Waliasing -Wampersand -Warray-bounds  \
           -Wcharacter-truncation -Wconversion -Wline-truncation -Wintrinsics-std  \
           -Wsurprising -Wno-tabs -Wunderflow -Wintrinsic-shadow \
           -Wno-align-commons -Werror -Wunused-parameter -Wextra
LDFLAGS  =  -fPIC
########################################################################################################
# Libraries to Link
LIB     = -ljmfft_i8r8
LIBDIR  = -L/home/fdubuffe/lib/
########################################################################################################
#
SRC      = mod_param.f90 inv_mat_band.f90 mod_figure.f90 mod_bounds.f90 \
           mod_adi.f90 mod_adi.f90 mod_NS.f90 mod_initio.f90 NS_code2D.f90
OBJ      = $(SRC:.f90=.o)

main: $(OBJ)
	$(F90) -o $@ $^ $(LDFLAGS) $(LIBDIR) $(LIB)

%.o: %.f90
	$(F90) -o $@ -c $< $(F90FLAGS)

cleanall: clean
	rm -f *_gmt gmt_exe *.grd *.cpt *.tmoy fort* *.ps *.epsi tp.gif tp*.plt

clean:
	rm -f *.o *.mod main fort.*
