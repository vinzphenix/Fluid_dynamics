
INC_DIR := -I/home/vincelinux/Documents/Fluid_dynamics/Project/lib_petsc/include
LIB_DIR := -L/home/vincelinux/Documents/Fluid_dynamics/Project/lib_petsc/lib -Wl,-rpath=/home/vincelinux/Documents/Fluid_dynamics/Project/lib_petsc/lib
LIB := -lpetsc

#INC_DIR := -I/usr/local/opt/petsc/include
#LIB_DIR := -L/usr/local/opt/petsc/lib
#LIB := -lpetsc.3.14

CXX_FLAGS := -J -O3 -Wall #-Werror #-g

#Compilation
all :
	gcc -o project project.c poisson.c -lm -lmpi $(CXX_FLAGS) $(LIB_DIR) $(LIB) $(INC_DIR)

#Delete of executable files
clean :
	rm projet

#Delete of results
clean_txt :
#	rm -vf results/*.txt results/P-* results/U-* results/V-* results/Reh-* results/Vtx-* results/Rehw-* results/Div-*
	rm data/*.txt
