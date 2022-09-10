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
The Wigner-Seitz cell is defined as the set of points which are closer to the lattice point $00$ than they are to any other lattice point. This set is a primitive unit cell because it contains just 1 lattice point ($00$) and can completely tile the lattice plane with no gaps nor overlaps. There is a very simple construction which can be used to find the vertices of this cell, which is explained below.

To tackle this problem, we first want to break it up into smaller manageable pieces. So, suppose we have two points: $00$, and some other point $P$. How can we find all points which are closer to the origin than they are to $P$? The answer is as follows:
1. Draw the line between the origin and $P$.
2. Draw the perpendicular bisector of this line.
3. All of the points which are on the origin's side of the line are closer to the origin, while all points on $P$'s side of the line are closer to $P$. This might seem obvious, but it will be very helpful when we scale up our problem.

Suppose now that we have three points: the origin, some point $P$, and another point $Q$. If we want to find all the points that are closer to the origin than they are to either $P$ *or* $Q$, we will basically do the same procedure we did above. We'll first find the set of points which are closer to the origin than $P$, and then we'll find the set of points that are closer to the origin than $Q$. Then, we'll take the intersection of both sets to get the points that are closer to the origin than either $P$ or $Q$.

Our lattice has an infinite number of sets, but we want a finite number of sets so we can carry out our procedure. Therefore, we need to find some set of lattice points that are "representative" of all of the lattice points, in a certain manner of speaking.

To do this, we will notice that if we have the point $P = m_1a + n_1b$, where $a$ and $b$ are lattice vectors and $m_1$ and $n_1$ are scalars, we do not also need any other point $Q = m_2a + n_2b$, where $m_2 > m_1$ and $n_2 > n_1$, because all the points closer to the origin than $P$ are also closer to the origin than $Q$, but some of the points closer to the origin than $Q$ are not closer to the origin than $P$. Therefore, $P$ is "more useful" in some sense for our purposes.

From this, we can see that we just need to check the lattice points that are closest to the origin. These eight points are:
- $\overline{1}\overline{1}$
- $0\overline{1}$
- $1\overline{1}$
- $10$
- $11$
- $01$
- $\overline{1}1$
- $\overline{1}0$

If we find the sets of points that are closer to the origin than to each of these points respectively, and then take the intersection, we get our desired set of points.

To get the vertices of this polygon, we'll just pairwise intersect the bisectors each of our points generated.

The user can generate and plot the Wigner-Seitz cell for any lattice, as well as find its area.

## Future Features
- Custom motif at each lattice point
- Lattice transformations
- Lattice symmetries
- Loading/saving views
- Area of standard primitive lattice cells