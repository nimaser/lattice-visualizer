# @author nimaser
# @version 0.1.0
# @date Sep 9, 2022

import itertools

import tkinter as tk
import matplotlib as mpl
import numpy as np

from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)

def polygon_area(x : list[float], y : list[float]) -> float:
    """
    Given the vertices of a polygon (defined in two lists of `x` coords and `y` coords), finds the
    area of a polygon using the Shoelace formula. From https://stackoverflow.com/a/30408825
    """
    # TODO: understand how this works
    return 0.5*np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))

def get_int_points_in_polygon(vertices : np.ndarray) -> np.ndarray:
    """
    Given a list of vertices defining a polygon, returns all points with integer coordinates within it.
    """
    # decompose the points into lists of x and y values
    x_vals = [vertex[0] for vertex in vertices]
    y_vals = [vertex[1] for vertex in vertices]

    # find the minimum and maximum x and y values
    x_min, x_max = min(x_vals), max(x_vals)
    y_min, y_max = min(y_vals), max(y_vals)

    # round the min and max values to the nearest int, adding 1 to the max for use in range()
    x_range = map(int, (np.ceil(x_min), np.floor(x_max + 1)))
    y_range = map(int, (np.ceil(y_min), np.floor(y_max + 1)))

    # get all possible points with integer coordinates inside of the rectangular region which most
    # tightly bounds the provided polygon by finding the cartesian product of the x and y ranges
    candidates = itertools.product(range(*x_range), range(*y_range))

    # a mpl object which can be used to find whether points are inside of a polygon
    polygon = mpl.path.Path(vertices)

    # boolean array determining which candidates are inside of the polygon
    is_inside = polygon.contains_points(candidates)

    # get a list of the candidates which are inside of the polygon via array compression
    return [candidates[i] for i in range(len(candidates)) if is_inside[i]]


class Lattice:
    def __init__(self, dimension : int, basis : np.ndarray):
        """
        The `Lattice` class represents a lattice with dimension `dimension` defined by `basis`, a set of column basis vectors.

        Points in space can be represented in two ways, either as a coordinate in standard space, which is defined by
        the standard unit vectors, or as a coordinate in lattice space, which is defined by `basis`. This class provides
        methods to convert between the two representations.

        Points whose lattice coordinate representation consists of only integers are considered valid lattice points.
        """
        assert len(basis) == dimension, f"Number of basis vectors must match dimension. {len(basis)} vectors with dimension {dimension}."
        for vector in basis:
            assert(len(vector)) == dimension, f"Length of basis vector must match dimension. {vector} has {len(vector)} coordinates with dimension {dimension}."

        self.dimension = dimension
        self.basis = np.array(basis)
        self.mat_lat_to_std = self.basis.transpose()
        self.mat_std_to_lat = np.linalg.inv(self.basis.transpose())
        
    def get_dimension(self):
        """Returns the dimension of this lattice."""
        return self.dimension

    def get_basis(self):
        """Returns the basis vectors of this lattice."""
        return self.basis

    def trans_lat_to_std(self, lat_points : list[np.ndarray]):
        """Transforms a point in lattice coordinates to its standard coordinates representation."""
        return [self.mat_lat_to_std @ point for point in lat_points]

    def trans_std_to_lat(self, std_points : list[np.ndarray]):
        """Transforms a point in standard coordinates to its lattice coordinates representation."""
        return [self.mat_std_to_lat @ point for point in std_points]

    def valid_lattice_points(self, std_points):
        """
        Returns whether the lattice representation of `std_point`, given in standard coordinates,
        consists of only integers. Due to the transformation required, some cases may fail due to
        floating-point issues.
        """
        # if a number equals its floor, it's an integer
        lat_points = np.array(self.trans_std_to_lat(std_points))
        return [all(point) for point in lat_points == np.floor(lat_points)]
    
    def validate_basis(self, other_basis):
        raise NotImplementedError()
    

class Lattice_2D(Lattice):
    def __init__(self, basis):
        super().__init__(2, basis)


class Lattice_3D(Lattice):
    def __init__(self, basis):
        super().__init__(3, basis)

def main():
    # initialize window
    WIDTH, HEIGHT = 1280, 720
    window = tk.Tk()
    window.title("lattice-visualizer")
    window.geometry(f"{WIDTH}x{HEIGHT}")
    window.resizable(False, False)

    # set up window's grid layout
    window.rowconfigure(0, weight=1)
    window.columnconfigure(0, weight=3)
    window.columnconfigure(1, weight=22)
    window.columnconfigure(2, weight=75)

    # set up and grid the three main containers
    sidebar = tk.Frame(window)
    sidebar.grid_propagate(False)
    sidebar.grid(row=0, column=0, sticky=tk.NSEW)

    config = tk.Frame(window, bg="blue")
    config.grid_propagate(False)
    config.grid(row=0, column=1, sticky=tk.NSEW)

    graph = tk.Frame(window)
    graph.grid_propagate(False)
    graph.grid(row=0, column=2, sticky=tk.NSEW)

    # set up the sidebar
    num_lattices = 4
    for i in range(num_lattices): sidebar.rowconfigure(i, weight=1)
    sidebar.columnconfigure(0, weight=1)

    lattice_buttons = []
    for i in range(num_lattices):
        lattice_buttons.append(tk.Button(sidebar, text=f"{i + 1}"))
        lattice_buttons[i].grid(row=i, column=0, sticky=tk.NSEW)

    # set up the config
    # title: Lattice #

    # set up the graph
    graph.rowconfigure(0, weight=1)
    graph.columnconfigure(0, weight=1)

    fig = mpl.figure.Figure(dpi=100)
    fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    ax = fig.add_subplot(111)
    ax.grid(linestyle="dashed")

    figcanvas = FigureCanvasTkAgg(fig, master=graph)
    figcanvas.draw()
    figcanvas.get_tk_widget().grid(row=0, column=0, sticky=tk.NSEW)

    # toolbar = NavigationToolbar2Tk(fig, window)
    # toolbar.update()
    # toolbar.get_tk_widget().pack()
    
    window.mainloop()


if __name__ == "__main__":
    main()