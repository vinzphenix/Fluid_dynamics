
INC_DIR := -Ilib_petsc/include
LIB_DIR := -L/home/vincelinux/Documents/Fluid_dynamics/Project/lib_petsc/lib -Wl,-rpath=/home/vincelinux/Documents/Fluid_dynamics/Project/lib_petsc/lib
LIB := -lpetsc

SRC_DIR := src
OBJ_DIR := obj
BIN_DIR := .

EXE := $(BIN_DIR)/cfd

SRC := $(wildcard $(SRC_DIR)/*.c)
OBJ := $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC))

CPPFLAGS := -Iinclude -MMD -MP # -I is a preprocessor flag, not a compiler flag
CFLAGS   := -O3 -Wall          # some warnings about bad code
LDFLAGS  := -Llib_petsc/lib -Wl,-rpath=lib_petsc/lib        # -L is a linker flag
LDLIBS   := -lm -lmpi -lpetsc         # Left empty if no libs are needed


.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJ) | $(BIN_DIR)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

# $(EXE): $(OBJ) | $(BIN_DIR)
# 	$(CC) $^ $(LDLIBS) $(LIB_DIR) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(INC_DIR) -c $< -o $@

$(BIN_DIR) $(OBJ_DIR):
	mkdir -p $@

# clean:
# 	@$(RM) -rv $(BIN_DIR) $(OBJ_DIR)

clean:
	@$(RM) -rv $(OBJ_DIR)

-include $(OBJ:.o=.d)




# CXX_FLAGS := -O3 -Wall #-Werror #-g

# #Compilation
# all :
# 	gcc -o cfd main.c project.c poisson.c adi.c -lm -lmpi $(CXX_FLAGS) $(LIB_DIR) $(LIB) $(INC_DIR)

# run :
# 	gcc -o cfd main.c project.c poisson.c adi.c -lm -lmpi $(CXX_FLAGS) $(LIB_DIR) $(LIB) $(INC_DIR)
# 	./cfd -ksp_type fgmres -pc_type lu


# #Delete of executable files
# clean :
# 	rm cdf

# #Delete of results
# clean_txt :
# #	rm -vf results/*.txt results/P-* results/U-* results/V-* results/Reh-* results/Vtx-* results/Rehw-* results/Div-*
# 	rm results/*.txt
