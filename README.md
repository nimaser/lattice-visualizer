# lattice-visualizer by nimaser

## Overview
This application provides basic plotting and manipulation functionality for planar (2D) crystallographic lattices.

## Background
A familiarity with basic linear algebra is assumed, so that the only new material is the specifics of planar crystallographic lattices.

Consider two vector spaces, the standard space $S$ defined by the familiar basis vectors $e1 = (1, 0)$ and $e2 = (0, 1)$, and the lattice space $L$ defined by two linearly independent custom basis vectors $a$ and $b$.

Any point on the plane has two coordinate-pair representations, one in $S$ and one in $L$. We can define a transformation that maps one space onto the other, allowing us to convert points in the standard space to their representations in the lattice space, and vice versa. This transformation can be represented as a matrix containing the basis vectors of the space.

We define $M$ to be the matrix that converts from lattice coordinates to standard coordinates, with $M = [a, b]$, and its inverse $M^-1$ is the matrix that converts from standard coordinates to lattice coordinates.

We will categorize all points in lattice space as being either "valid lattice points" or not, based on whether their lattice coordinate representation consists of only integers. For example, the point (1, -3) is a valid lattice point, while the point (1.5, 2.3) is not. Note that valid lattice points may have non-integer standard coordinate representations, and that is fine.

## Plotting a Lattice
The user first defines a lattice L by providing the two basis vectors $a$ and $b$ and an offset location $o$. The offset is simply the location in standard space where the lattice point (00) will be located, and all other lattice points will be calculated based on it.

Next, the user defines a region R in the standard space in one of two ways: either by providing a list of vertices which describe a polygon, or by providing a radius and center point for a circular region.

Finally, the user chooses a color for the lattice points. At this point, the application can then plot all valid lattice points of L which lie within R. The user can toggle on/off a display of the lattice's basis vectors.

A total of four different lattices can be plotted at one time.

## Axis Settings
The user can toggle gridlines on/off and set the viewing region (axis bounds).

## Graph Options
The user can left-click on lattice points to toggle highlighting them and displaying their lattice coordinates. The user can right-click on lattice points to toggle a vector decomposition for the point, which shows how the basis vectors can be added to obtain that point.

## Wigner-Seitz Cell
The user can step through the steps to find the Wigner-Seitz cell for a lattice. The steps are explained below:
TODO

## Future Features
- Custom motif at each lattice point
- Lattice transformations
- Lattice symmetries
- Loading/saving views