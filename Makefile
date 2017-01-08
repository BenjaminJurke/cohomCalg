#########################################################################
#                                                                       #
#  Makefile for the cohomCalg application (including modified PolyLib)  #
#                                                                       #
#########################################################################

# compilers & flags
CC        := g++
CFLAGS    := -O3
LD        := g++
LDFLAGS   := 
# Note: (1) If you want to enforce 32-bit or 64-bit compilation, 
#           add "-m32" or "-m64" to both CFLAGS and LDFLAGS.
#       (2) If you want to enforce a static linking of all libraries 
#           (i.e. include everything into one binary), add "-static" to LDFLAGS


# global defs
DEFS      := -DPOLYLIB_BITS=64

# directories
COHOMCALG_SRC_DIR := source
POLYLIB_SRC_DIR   := source/polylib_mod
BUILD_DIR         := build/polylib_mod build

# source and object files
COHOMCALG_SRC := $(foreach sdir,$(COHOMCALG_SRC_DIR),$(wildcard $(sdir)/*.cpp))
COHOMCALG_OBJ := $(patsubst source/%.cpp,build/%.o,$(COHOMCALG_SRC))
POLYLIB_SRC   := $(foreach sdir,$(POLYLIB_SRC_DIR),$(wildcard $(sdir)/*.c))
POLYLIB_OBJ   := $(patsubst source/%.c,build/%.o,$(POLYLIB_SRC))
INCLUDES      := $(addprefix -I,$(COHOMCALG_SRC_DIR) $(POLYLIB_SRC_DIR))

vpath %.c $(POLYLIB_SRC_DIR)
vpath %.cpp $(COHOMCALG_SRC_DIR)

# macros for .c/.cpp dirs
define make-goal-cpp
$1/%.o: %.cpp
	$(CC) $(INCLUDES) $(DEFS) $(CFLAGS) -c $$< -o $$@
endef

define make-goal-c
$1/%.o: %.c
	$(CC) $(INCLUDES) $(DEFS) $(CFLAGS) -c $$< -o $$@
endef


.PHONY: all checkdirs clean

all: checkdirs bin/cohomcalg

bin/cohomcalg: $(COHOMCALG_OBJ) $(POLYLIB_OBJ)
	$(LD) $^ $(LDFLAGS) -lpthread -o $@
#-static -lsource/polylib-5.22.5/.libs/libpolylib64.a

checkdirs: $(BUILD_DIR)

$(BUILD_DIR):
	@mkdir -p $@

clean:
	@rm -rf $(BUILD_DIR)


$(eval $(call make-goal-c, build/polylib_mod))
$(eval $(call make-goal-cpp, build))