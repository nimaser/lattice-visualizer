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


class Lattice2D:
    def __init__(self, basis : np.ndarray, offset : np.ndarray):
        """
        Points in space can be represented in two ways, either as a coordinate in standard space,
        which is defined by the standard unit vectors, or as a coordinate in lattice space, which
        is defined by `basis`. This class provides methods to convert between the two.

        The `Lattice2D` class represents a 2D lattice defined by `basis`, a set of column basis
        vectors, and `offset`, the position in the standard space of the 00 lattice point.

        Points whose lattice coordinate representation consists of only integers are considered
        valid lattice points.
        """
        # Safety first!
        assert len(basis) == 2, f"{len(basis)} basis vectors provided, while 2 required."
        for vector in basis:
            assert(len(vector)) == 2, f"Length of basis vector {vector} must be 2. Was {len(vector)}"

        self.basis = np.array(basis)
        self.mat_lat_to_std = self.basis.transpose()
        self.mat_std_to_lat = np.linalg.inv(self.basis.transpose())

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

    def get_

    def generate_wigner_seitz_lines(u : tuple[float, float], v : tuple[float, float]) -> list[tuple[tuple[float, float], tuple[float, float]]]:
        """
        Given the two basis vectors u and v of the lattice, returns a list of line segments which
        border the Wigner Seitz cell.
        """
        P_N_to_S = lambda point: transform_point_N_to_S(u, v, point)

        # first get the 8 lattice points closest to 00, in clockwise order from -1 -1
        x, y = [-1, -1, -1, 0, 1, 1, 1, 0], [-1, 0, 1, 1, 1, 0, -1, -1]
        transformed_x, transformed_y = [], []

        # iterate through the transformation of the cartesian product of {-1, 0, 1} with itself
        for point in zip(x, y):
            x1, y1 = P_N_to_S(point)
            transformed_x.append(x1)
            transformed_y.append(y1)

        # since the Wigner Seitz cell is guaranteed to be a hexagon for 2D nonrectangular lattices,
        # we only need 6 lines to define it. in the case of a rectangular lattice, we only need 4, and
        # the other four have redundant information. however, in the case of the hexagon, two of the
        # lines do not contribute to the hexagon, making processing them to find the vertices (and hence
        # the area) difficult. therefore, we want to remove them. looking at many examples, we see that
        # the lines created based on the two lattice points furthest from the origin are the problematic
        # ones. therefore, we will remove these two lattice points from our list
        for _ in range(2):
            max_distance = 0
            max_index = 0
            for i in range(len(transformed_x)):
                x1, y1 = transformed_x[i], transformed_y[i]
                d = np.sqrt(x1 ** 2 + y1 ** 2)
                if d > max_distance:
                    max_distance = d
                    max_index = i
            transformed_x.pop(max_index)
            transformed_y.pop(max_index)

        # each line will be a 2-tuple containing two 2-tuples, one for the x points and one for the y
        lines = []

        # cutting each lattice point in 2 gives the bisecting point. adding (-b, a) and (b, -a) to this
        # point, where (a, b) is one of the lattice points, gives endpoints for a line guaranteed to be
        # longer than the edge of the Wigner Seitz cell which it borders
        for a, b in zip(transformed_x, transformed_y):
            P1 = (a / 2 - b, b / 2 + a)
            P2 = (a / 2 + b, b / 2 - a)
            line = ((P1[0], P2[0]), (P1[1], P2[1]))
            lines.append(line)

        return lines

    def get_wigner_seitz_vertices(lines : list[tuple[tuple[float, float], tuple[float, float]]]) -> list[tuple[float, float]]:
        """
        Given a list of lines describing a Wigner-Seitz cell, returns the vertices of the cell. Assumes
        that the lines are given in sequential order, because it finds the vertices by finding the
        intersections of the lines. Uses code from https://stackoverflow.com/a/3252222
        """
        def perp(a) :
            b = np.empty_like(a)
            b[0] = -a[1]
            b[1] = a[0]
            return b

        def seg_intersect(a1, a2, b1, b2) :
            da = a2-a1
            db = b2-b1
            dp = a1-b1
            dap = perp(da)
            denom = np.dot( dap, db)
            num = np.dot( dap, dp )
            return (num / denom.astype(float))*db + b1

        vertices = []

        # intersect each line with the previous one, so 0 with -1, 1 with 0, etc to n - 1 with n - 2,
        # where n is the number of lines
        for i in range(len(lines)):
            a = lines[i - 1]
            b = lines[i]
            a1 = np.array((a[0][0], a[1][0]))
            a2 = np.array((a[0][1], a[1][1]))
            b1 = np.array((b[0][0], b[1][0]))
            b2 = np.array((b[0][1], b[1][1]))
            vertices.append(seg_intersect(a1, a2, b1, b2))

        # remove duplicate points in the case of a rectangular WS cell, from https://stackoverflow.com/a/44628266
        lookup = set()
        #vertices = [vertex for vertex in vertices if vertex not in lookup and lookup.add(vertex) is None]

        # split into x and y lists
        vx, vy = [], []
        for vertex in vertices:
            vx.append(vertex[0])
            vy.append(vertex[1])
        return vx, vy
    
class LatticePlot:
    def __init__():
        pass



class AxesManager:
    def __init__():
        pass

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
    for i in range(num_lattices):
        sidebar.rowconfigure(i, weight=1)
    sidebar.columnconfigure(0, weight=1)

    lattice_buttons = []
    for i in range(num_lattices):
        lattice_buttons.append(tk.Button(sidebar, text=f"{i + 1}"))
        lattice_buttons[i].grid(row=i, column=0, sticky=tk.NSEW)

    # set up the config
    config.rowconfigure
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