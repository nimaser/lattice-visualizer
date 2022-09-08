# lattice-visualizer

## Overview
A lattice is a mathematical structure which is very useful in the field of crystallography. This application provides basic plotting and manipulation functionality for planar (2D) lattices.

## Background
Consider the standard basis vectors $e1$ and $e2$. We can represent any point on the plane as a linear combination of these two vectors, so we say that $e1$ and $e2$ "span" the plane.

As a reminder, a linear combination of two vectors takes the following form: $P = au + bv$, where $a$ and $b$ are scalars, $u$ and $v$ are basis vectors, and $P$ is the resulting point. To span the plane, we let $a$ and $b$ be all real numbers. Now, suppose that we limited the values of $a$ and $b$ to be only integers. Now, we have all points on the plane with integer coordinates. This is a square planar lattice: a periodic arrangement of discrete points on a plane such that each point lands on the vertex of a square.

Suppose instead that we chose different basis vectors, so $u$ and $v$ are no longer $e1$ and $e2$ but are instead some other custom vectors. We could end up with a lattice whose points have non-integer coordinates and are the vertices of a rectangle or a parallelogram.

## Standard vs Lattice Space
We are mainly concerned with two spaces: "Standard Space" and "Lattice Space". Standard space is our familiar Cartesian plane defined by the two basis vectors $e1$ and $e2$. Lattice space is defined by two custom vectors $a$ and $b$ (not the same as the scalars from earlier).