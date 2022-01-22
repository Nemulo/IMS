# TUMOR GROWTH SIMULATION

## IMS project - cellular automaton

### Author: Marek Nemeth

## Build

use
```
make
```
to build program

## Running program

program accepts following options
```
./tumor -h
```
shows help

```
./tumor -p <probability> -R <max Radius> -a <factor a> -b <factor b> -g
```
shows GUI for user in window using GL library

### Requred arguments
 
```
-p <probability>
```
sets base probability for proliferate cells to divide - use float value between 0-1

```
-R <max Radius>
```
defines maximal radius of tumor - use float value, tumor radius is set in mm

```
-a <factor a>
```
Base  necrotic  thickness,  controlled  by  nutritional needs, use float value (there is required specific ratio between a factor and b factor to have proper results), set in mm

``` 
-b <factor b>
```
Base  proliferative  thickness,  controlled  by  nutritional needs, use float value (there is required specific ratio between a factor and b factor to have proper results), set in mm

## Required libraries

```
GL/glut.h
```

## Included files

in src/ directory
```
cell.cpp
cell.hpp
grid.cpp
grid.hpp
main.cpp
Makefile
```

documentation
```
README.md
documentation.pdf
```

