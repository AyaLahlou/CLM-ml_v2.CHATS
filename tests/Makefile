# ###############################################
# Makefile for CLM-ml functional test suite
#
# Prerequisites:
#   1. Build the main model first:
#        cd ../offline_executable && make
#   2. Then build test executables:
#        make all
#   3. Then run tests:
#        cd python && python -m pytest -v
#
# Environment variables required (same as main build):
#   NETCDF_PATH  - path to NetCDF installation
# ###############################################

MAIN_BUILD = ../offline_executable
OBJ        = $(MAIN_BUILD)/obj
LIB_NETCDF = $(NETCDF_PATH)/lib
MOD_NETCDF = $(NETCDF_PATH)/include

# Compiler: same flags as the main build
cmplr = nvfortran -C -Minform=inform -Ktrap=fp \
        -I$(MAIN_BUILD) -I$(MOD_NETCDF) \
        -L$(LIB_NETCDF) -lnetcdff -lblas -lm -L/usr/lib64

BIN = bin
SRC = fortran

# All main-build objects EXCEPT the main program entry point
MAIN_OBJS = $(filter-out $(OBJ)/CLMml.o, $(wildcard $(OBJ)/*.o))

# List of test executables
TARGETS = \
  $(BIN)/test_MLMathToolsMod_quadratic.exe              \
  $(BIN)/test_MLMathToolsMod_tridiag.exe                \
  $(BIN)/test_MLMathToolsMod_log_gamma_function.exe     \
  $(BIN)/test_MLMathToolsMod_beta_function.exe          \
  $(BIN)/test_MLMathToolsMod_beta_distribution_pdf.exe  \
  $(BIN)/test_MLMathToolsMod_beta_distribution_cdf.exe  \
  $(BIN)/test_MLWaterVaporMod_SatVap.exe                \
  $(BIN)/test_MLWaterVaporMod_LatVap.exe                \
  $(BIN)/test_MLLeafPhotosynthesisMod_ft.exe            \
  $(BIN)/test_MLLeafPhotosynthesisMod_fth.exe           \
  $(BIN)/test_MLLeafPhotosynthesisMod_fth25.exe         \
  $(BIN)/test_MLCanopyTurbulenceMod_phim_monin_obukhov.exe \
  $(BIN)/test_MLCanopyTurbulenceMod_psim_monin_obukhov.exe \
  $(BIN)/test_shr_orb_mod_shr_orb_params.exe            \
  $(BIN)/test_shr_orb_mod_shr_orb_decl.exe              \
  $(BIN)/test_MLMathToolsMod_tridiag_2eq.exe            \
  $(BIN)/test_MLLeafPhotosynthesisMod_RealizedRate.exe  \
  $(BIN)/test_MLCanopyTurbulenceMod_GetBeta.exe         \
  $(BIN)/test_MLCanopyTurbulenceMod_GetPrSc.exe

.PHONY: all clean

all: $(BIN) $(TARGETS)

$(BIN):
	mkdir -p $(BIN)

# Pattern rule: compile test driver and link against all main-build objects
$(BIN)/%.exe: $(SRC)/%.F90 $(MAIN_OBJS)
	$(cmplr) -o $@ $< $(MAIN_OBJS)

clean:
	rm -f $(BIN)/*.exe
