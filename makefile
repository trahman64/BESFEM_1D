# MFEM configuration
MFEM_INSTALL_DIR = /opt/software-current/2023.06/x86_64/generic/software/MFEM/4.7-foss-2023a

CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk

-include $(CONFIG_MK)

# Compiler and flags

CXX = $(MFEM_CXX)

CXXFLAGS = -g -O3 -std=c++17

INCLUDE_FLAGS = $(MFEM_INCFLAGS) -Iinclude
LIB_FLAGS     = $(MFEM_LIBS)


# Source files
SRC_FILES = src/battery_simulation.cpp src/linear_diffusion.cpp src/spherical_diffusion.cpp
    
#     src/diffusion.cpp \
#     src/radial_diffusion.cpp \
#     src/potential.cpp \    

# Executable

EXEC_DIR = bin
EXEC_NAME = solve1d
EXEC = $(EXEC_DIR)/$(EXEC_NAME)


# ====================================
# Default target
# ====================================

all: $(EXEC)


# Build

$(EXEC): $(SRC_FILES)
	@mkdir -p $(EXEC_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $^ -o $@ $(LIB_FLAGS)


# Clean

clean:
	rm -rf $(EXEC_DIR)

.PHONY: all clean