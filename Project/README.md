# 2D flow simulation - Project LMECA2660

This project was done in the course _Numerical methods in fluid mechanics_. It simulates a 2D flow around a rectangular obstacle. Many situations can be considered: the obstacle can move horizontally and/or vertically, the temperature field can be coupled with the velocity fields through natural convection under the Boussinesq approximation.


## Structure

    ├── README.md
    ├── Makefile
    │
    ├── anim              <- Animations in .mp4 format of different situations
    ├── doc               <- Documents, instructions, reference book
    ├── figures           <- Figures analyzing different simulations
    ├── include           <- Header files
    │   ├── project.h          # Physical parameters and solver options
    │   └── ...
    ├── lib_petsc         <- PETSc dependency
    ├── report            <- Report in .tex and .pdf formats
    ├── results           <- Results of simulations saved in .txt files
    ├── scripts           <- Post-processing of the simulations: animations and figures
    └── src               <- Source .c files

## Run the code
The code in __C__ can be compiled easily with make, and then executed as follows:
```
make -j
./cfd -ksp_type fgmres -pc_type lu -n 50 -dt 0.002 -tend 50. -freq 0.1 -dir new_case
```

- __-ksp_type__ : Name of the PETSc Krylov method, that shouldn't be modified
- __-pc_type__ : Name of the PETSc preconditioner method, that shouldn't be modified
- __-n__ [integer] : Spatial discretization, with $H_{box}=1/n$
- __-dt__ [double] : Time step, or initial time step if adaptative time step activated
- __-tend__ [double] : Final time of the simulation
- __-freq__ [double] : Frequency to which the program saves the fields in `.txt` files
- __-dir__ : Name of the subdirectory created or overwritten in `./results/` directory


## Produce an animation
The animation proposed shows the streamlines (iso-$\psi$), the pressure field $p$, the vorticity field $\omega$ and if possible the temperature field $T$.

First, you need to produces the frames of the animation based on a simulation stored in the directory `./results/`. To do so, make sure that the directory `./anim/` exists and execute:
```
python3 run.py -dir new_case
```

Then, to produce a `.mp4` file, you can use `FFmpeg` as follows:
```
ffmpeg -framerate 24 -i ./anim/new_case/frame_%05d.png -vcodec libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -crf 5 -r 25 -pix_fmt yuv420p ./anim/new_case.mp4
```

## Dependencies
The solver uses the PETSc library in order to solve the poisson equation of the reprojection scheme. Instructions for installing are available in the `./doc/` directory, or on the website https://petsc.org/release/install/install/.

## Schema of the compuatational domain with default dimensions.
```
────────────────────────────────────────────────────────────────
->                        5                                     
->              <------------------->                           
->                                                              
->      3      ╭─────────────────────╮             7            
->  <------->  │              | 1    │  <---------------------> 
->             ╰─────────────────────╯                          
->                        |                                     
->                        | 2                                   
->                        |                                     
────────────────────────────────────────────────────────────────
```