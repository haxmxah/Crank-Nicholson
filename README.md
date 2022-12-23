# Crank-Nicholson
Resolution of EDP, parabolic equations, Poisson equation 1D and heat equation via Crank Nicholson.

## Diffusion process
Governed by the equation

$$\frac{\partial \theta(x,t)}{\partial t}= \alpha \frac{\partial^2 \theta (x,t)}{\partial x^2}-\beta \theta(x,t)$$

We study its behavior.

# Compilation and execution of the program
This program was written in _Fortran_ 77 and the graphics were plotted with _Gnuplot_.
## Linux and Mac
### Compilation

```
gfortran -name_of_the_file.f -o name_of_the_output_file.out
```
### Execution
```
./name_of_the_output_file.out
```

## Windows
### Compilation
```
gfortran -name_of_the_file.f -o name_of_the_output_file.exe
```
### Execution
```
./name_of_the_output_file.exe
```
