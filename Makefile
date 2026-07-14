FC      = gfortran

FFLAGS = -Wall -Wextra -fcheck=all -g \
          -fopenmp \
          -Jbuild/mod \
          -Ibuild/mod

# =====================================================
# HEALPix
# =====================================================

HEALPIX = /home/pbuser/Schreibtisch/healpix_fortran/Healpix_3.83

INCLUDES = -I$(HEALPIX)/include
LIBDIR   = -L$(HEALPIX)/lib

LIBS = \
    -lblas \
    -llapack \
    -lhealpix \
    -lsharp \
    -lcfitsio

# =====================================================

SRC_DIR = build/src
MOD_DIR = build/mod
OBJ_DIR = build/obj

TARGET  = run/mdgraph

# =====================================================
# SOURCE FILES
# =====================================================

SRC := \
    $(SRC_DIR)/core/topology_module.f90 \
    $(SRC_DIR)/io/dcd_header_reader.f90 \
    $(SRC_DIR)/core/project_module.f90 \
    $(SRC_DIR)/core/analysis_plan_module.f90 \
    $(SRC_DIR)/core/time_window_file_mapping.f90 \
    $(SRC_DIR)/core/test_time_window_file_mapping.f90 \
    $(SRC_DIR)/io/adaptive_dcd_reader.f90 \
    $(SRC_DIR)/io/io_pdb_module.f90 \
    $(SRC_DIR)/io/pdb_read.f90 \
    $(SRC_DIR)/io/dcd_reader.f90 \
    $(SRC_DIR)/core/time_window_loader.f90 \
    $(SRC_DIR)/core/adaptive_stream_module.f90 \
    $(SRC_DIR)/features/orientation_vectors.f90 \
    $(SRC_DIR)/core/analysis_executor_module.f90 \
    $(SRC_DIR)/core/system_module.f90 \
    $(SRC_DIR)/core/sidechains.f90 \
    $(SRC_DIR)/core/domains_module.f90 \
    $(SRC_DIR)/core/domain_refstruc.f90 \
    $(SRC_DIR)/core/sidechain_tensors.f90 \
    $(SRC_DIR)/core/tra_lift_module.f90 \
    $(SRC_DIR)/core/sidechainhead_refstruc.f90 \
    $(SRC_DIR)/core/direction_vectors_analyse_plan_module.f90 \
    $(SRC_DIR)/core/position_lift_module.f90 \
    $(SRC_DIR)/core/healpix_type_module.f90 \
    $(SRC_DIR)/core/healpix_sampling_plan_module.f90 \
    $(SRC_DIR)/core/healpix_marginals.f90 \
    $(SRC_DIR)/core/healpix_transfer_entropy_module.f90 \
    $(SRC_DIR)/core/healpix_analysis_execution_module_single_lag.f90 \
    $(SRC_DIR)/core/main.f90

# =====================================================
# OBJECT FILES
# =====================================================

OBJ := $(SRC:$(SRC_DIR)/%.f90=$(OBJ_DIR)/%.o)

# =====================================================
# DEFAULT
# =====================================================

all: dirs $(TARGET)

# =====================================================
# DIRECTORIES
# =====================================================

dirs:
	@mkdir -p $(MOD_DIR)
	@mkdir -p $(OBJ_DIR)/core
	@mkdir -p $(OBJ_DIR)/io
	@mkdir -p $(OBJ_DIR)/features
	@mkdir -p run

# =====================================================
# FULL CLEAN BUILD
# =====================================================

clean_build:
	rm -rf $(OBJ_DIR)
	rm -rf $(MOD_DIR)

# =====================================================
# LINKING
# =====================================================

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o $(TARGET) $(LIBDIR) $(LIBS)
# =====================================================
# COMPILATION RULE
# =====================================================

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) $(INCLUDES) -c $< -o $@

# =====================================================
# RUN
# =====================================================

run: all
	./$(TARGET)

# =====================================================
# CLEAN
# =====================================================

clean:
	rm -rf $(OBJ_DIR) $(MOD_DIR)
	rm -f $(TARGET)

# =====================================================
# DEBUG
# =====================================================

list:
	@echo "Source files:"
	@for f in $(SRC); do echo $$f; done