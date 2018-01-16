# Projects in Computational Physics
This repository stores a number of projects completed as part of a course in
Computational Physics, studied at the University of Melbourne in 2016 as the
course PHYC30012.

This repository stores five projects. Each required calculations or simulations
to be written in C, along with an extended written report.

### Project 1 - Semiclassical Quantisation of Molecular Vibrations
Calculation of the allowed energy levels of vibration for molecular hydrogen,
using these to then investigate the Morse potential as a model for the actual
potential of molecular hydrogen.

The associated program first compares the trapezoidal and Simpson methods for
numerical integral evaluation. It then uses Simpson's method and the false
position algorithm to determine the root of an integral equation. Putting these
functions into practice, it then determines the molecular vibrational energy
levels for two different potentials. First, a quadratic potential that can be
evaluated analytically, to test the program and gauge its accuracy. Then, the
Morse potential, a far more realistic function. The program outputs the
calculated energy levels for both potentials, and the known associated energies.
It optimises the variable *a* in the Morse potential to fit the first and the
second energy levels of molecular hydrogen, and contrasts how the other energy
levels drift as *a* is varied.

The program can be compiled by executing:

```
gcc -o project1 project1.c
```

It can by run by executing:

```
./project1
```

This uses default values for the number of steps to use, starting point,
stepsize, and tolerance. These can instead be entered on the command line, for
example:

```
./project1 100000 .1 .1 1.0e-7
```

![Project 1 expected output](assets/gifs/project1.gif)
