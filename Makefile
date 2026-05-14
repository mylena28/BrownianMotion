# =============================================================================
# Makefile – Brownian Dynamics: Settling Sphere
# Usage:
#   make            -> build both executables
#   make run        -> build + run the settling simulation
#   make validate   -> build + run the MSD validation
#   make all_run    -> build + run both
#   make clean      -> remove all build artefacts and output files
# =============================================================================

FC      = gfortran
FFLAGS  = -O2 -Wall -Wextra -fcheck=bounds

# Module source files (must be compiled before the programs that USE them)
MOD_SRCS = params_mod.f90 rng_mod.f90
MOD_OBJS = $(MOD_SRCS:.f90=.o)

# --- Default target ----------------------------------------------------------
.PHONY: all run validate all_run clean

all: brownian_sim msd_validate

# --- Executables -------------------------------------------------------------
brownian_sim: $(MOD_OBJS) brownian_sim.o
	$(FC) $(FFLAGS) -o $@ $^

msd_validate: $(MOD_OBJS) msd_validate.o
	$(FC) $(FFLAGS) -o $@ $^

# --- Generic rule: .f90 -> .o -----------------------------------------------
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# --- Inter-module dependencies (explicit) ------------------------------------
brownian_sim.o: params_mod.o rng_mod.o
msd_validate.o: params_mod.o rng_mod.o
rng_mod.o:      params_mod.o

# --- Convenience run targets -------------------------------------------------
run: brownian_sim
	./brownian_sim

validate: msd_validate
	./msd_validate

all_run: brownian_sim msd_validate
	./brownian_sim
	./msd_validate

# --- Clean -------------------------------------------------------------------
clean:
	rm -f *.o *.mod brownian_sim msd_validate
	rm -f trajectory.dat no_gravity_traj.dat msd_output.dat
	rm -f *.png *.mp4
