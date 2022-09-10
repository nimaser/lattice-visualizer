# @author nimaser
# @version 0.1.0
# @date Sep 9, 2022

import itertools

import tkinter as tk
import matplotlib as mpl
import numpy as np

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.animation import FuncAnimation

NUM_LATTICES = 4

def polygon_area(vertices : np.ndarray) -> float:
        """
        Given the vertices of a polygon, finds the area of a polygon using the Shoelace formula.
        From https://stackoverflow.com/a/30408825
        """
        x, y = zip(*vertices) # interesting fact: zip is basically its own inverse
        return 0.5*np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))

def get_circle_bbox(center, radius):
    """Returns the vertices of a rectangle that most closely contains the specified circle."""
    x_min = center[0] - radius
    x_max = center[0] + radius
    y_min = center[1] - radius
    y_max = center[1] + radius
    return (x_min, y_min), (x_max, y_min), (x_max, y_max), (x_min, y_max)

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
    candidates = list(itertools.product(range(*x_range), range(*y_range)))

    # a mpl object which can be used to find whether points are inside of a polygon
    polygon = mpl.path.Path(vertices)

    # boolean array determining which candidates are inside of the polygon
    is_inside = polygon.contains_points(candidates)

    # get a list of the candidates which are inside of the polygon via array compression
    return [candidates[i] for i in range(len(candidates)) if is_inside[i]]
    
class LatticePlot:
    
    # bounds enum
    POLYGONAL_BOUNDS = 0
    RADIAL_BOUNDS = 1

    # zordering:
    # 0: Wigner-Seitz cell
    # 1: Wigner-Seitz lines
    # 2: Wigner-Seitz intersection points
    # 3: Basis vectors
    # 4: Lattice points
    # 5: Axis Origin
    
    def __init__(self, ax):
        """Initializes the lattice to all None values."""
        self.ax = ax
        self.first_time_setup = True

        self.basis = None
        self.offset = None
        self.mat_lat_to_std = None
        self.mat_std_to_lat = None
        self.bounds = None
        self.bounds_type = None
        self.point_color = None

    def generate_lattice(self):
        # Calculate lattice points
        if self.bounds_type == LatticePlot.RADIAL_BOUNDS:
            vertices = get_circle_bbox(self.bounds["center"], self.bounds["radius"])
            vertices = self.trans_std_to_lat(vertices)
            candidates = get_int_points_in_polygon(vertices)
            lattice_points = []
            for point in candidates:
                if np.linalg.norm(np.array(point) - self.offset) < self.bounds["radius"]:
                    lattice_points.append(point)
            lattice_points = self.trans_lat_to_std(lattice_points)
        elif self.bounds_type == LatticePlot.POLYGONAL_BOUNDS:
            vertices = self.trans_std_to_lat(self.bounds)
            lattice_points = get_int_points_in_polygon(vertices)
            lattice_points = self.trans_lat_to_std(lattice_points)
        else: raise Exception("No bounds type selected.")
        self.lattice_points = lattice_points

        self.lattice = self.ax.plot(*zip(*lattice_points), f"{self.color} .")

        # if self.first_time_setup:
        #     #v1 = self.ax.quiver(0, 0, *self.basis[0], color='turquoise', scale_units='xy', scale=1, zorder=3)
        #     #v2 = self.ax.quiver(0, 0, *self.basis[1], color='orange', scale_units='xy', scale=1, zorder=3)
        #     #self.basis_arrows = [v1, v2]
        #     self.lattice_points = self.ax.plot(, zorder=4)[0]
            
        #     #self.ws_points = self.ax.plot([np.nan], [np.nan], zorder=2)[0]
        #     self.first_time_setup = False
        # else:
        #     # create basis vectors
        #     #self.basis_arrows[0].set_UVC(*self.basis[0], "b")
        #     #self.basis_arrows[1].set_UVC(*self.basis[1], "r")

        #     # create lattice points

        #     # create Wigner-Seitz cell
        #     pass

    def set_basis(self, basis):
        """Sets the basis of this lattice."""
        self.basis = np.array(basis)
        self.mat_lat_to_std = self.basis.transpose()
        self.mat_std_to_lat = np.linalg.inv(self.basis.transpose())

    def set_offset(self, offset):
        """Sets the 00 offset of this lattice."""
        self.offset = np.array(offset)

    def set_color(self, color):
        """Sets the color of this lattice."""
        self.color = color

    def set_bounds_polygonal(self, vertices):
        """Sets the bounds of this plotted lattice using a polygonal scheme."""
        self.bounds = vertices
        self.bounds_type = LatticePlot.POLYGONAL_BOUNDS

    def set_bounds_radial(self, radius, center):
        """Sets the bounds of this plotted lattice using a radial scheme."""
        self.bounds = {"radius" : radius, "center" : center}
        self.bounds_type = LatticePlot.RADIAL_BOUNDS

    def set_point_color(self, color):
        """Sets the color of this lattice's points."""
        self.color = color

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
        # if n = floor(n), n is an int
        lat_points = np.array(self.trans_std_to_lat(std_points))
        return [all(point) for point in lat_points == np.floor(lat_points)]

    def get_unit_cell_area(self):
        """Returns the area of the primitive unit cell."""
        return np.linalg.norm(np.cross(*self.basis))

    def update(self):
        pass

    def generate_wigner_seitz_lines(self):
        """
        Returns a list of lines which border the Wigner Seitz cell.
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


class AxesManager:
    def __init__(self, fig, toolbar):
        # Configure toolbar
        self.toolbar = toolbar
        self.toolbar._Spacer()
        self.toolbar._Button("Toggle grid lines", None, True, self.toggle_grid_visibility)
        self.toolbar._Spacer()
        self.toolbar._Button("Toggle origin", None, True, self.toggle_origin)
        # Set up figure
        self.fig = fig
        self.fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
        # Set up axes, grid, and origin point
        self.ax = fig.add_subplot(111)
        self.ax.grid(linestyle="dashed")
        self.grid_visibility = False
        self.ax.grid(visible=self.grid_visibility)
        self.origin_visibility = False
        self.origin_point = self.ax.plot([0], [0], "k .", zorder=5, visible=self.origin_visibility)[0]
        # Set up animation
        self.ani = FuncAnimation(fig, self.update, blit=False)

        # Set up lattices
        self.lattice_plots = []
        for i in range(NUM_LATTICES):
            self.lattice_plots.append(LatticePlot(self.ax))
        
        # testing zone
        self.lattice_plots[0].set_basis([(1, 0), (0, 1)])
        self.lattice_plots[0].set_offset([0, 0])
        self.lattice_plots[0].set_color("b")
        self.lattice_plots[0].set_bounds_radial(3, (0, 0))
        self.lattice_plots[0].generate_lattice()

        self.lattice_plots[1].set_basis([(1, 1), (-1, 1)])
        self.lattice_plots[1].set_offset([0, 0])
        self.lattice_plots[1].set_color("r")
        self.lattice_plots[1].set_bounds_radial(5, (0, 0))
        self.lattice_plots[1].generate_lattice()

    def toggle_grid_visibility(self):
        self.grid_visibility = not self.grid_visibility
        self.ax.grid(visible=self.grid_visibility)

    def toggle_origin(self):
        self.origin_visibility = not self.origin_visibility
        self.origin_point.set_visible(self.origin_visibility)

    def update(self, frame):
        pass


class LatticeConfig():
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

    # set up the graph, figure, toolmbar, and manager
    graph.rowconfigure(0, weight=1)
    graph.columnconfigure(0, weight=1)
    fig = mpl.figure.Figure(dpi=100)
    figcanvas = FigureCanvasTkAgg(fig, master=graph)
    figcanvas.draw()
    figcanvas.get_tk_widget().grid(row=0, column=0, sticky=tk.NSEW)
    toolbar_frame = tk.Frame(graph)
    toolbar_frame.grid(row=1, column=0, sticky=tk.NSEW)
    toolbar = NavigationToolbar2Tk(figcanvas, toolbar_frame)
    axes_manager = AxesManager(fig, toolbar)

    # set up the config
    # config.rowconfigure
    # title: Lattice #

    # set up the sidebar
    for i in range(NUM_LATTICES):
        sidebar.rowconfigure(i, weight=1)
    sidebar.columnconfigure(0, weight=1)

    lattice_buttons = []
    for i in range(NUM_LATTICES):
        lattice_buttons.append(tk.Button(sidebar, text=f"{i + 1}"))
        lattice_buttons[i].grid(row=i, column=0, sticky=tk.NSEW)

    window.mainloop()


if __name__ == "__main__":
    main()