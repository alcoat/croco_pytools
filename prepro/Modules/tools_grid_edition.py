#--- Dependencies ---------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import distance
from ipywidgets import widgets, FloatText, VBox, Button, interact_manual, BoundedIntText, interactive, IntSlider
from IPython.display import display
from matplotlib.widgets import RectangleSelector

#--- Functions / Widgets ---------------------------------------------------

def plot_outline_ax(outputs, ax):
    def outline(lon, lat):
        '''
        Return lon, lat of perimeter around the grid
        '''
        def func(var):
            return np.hstack([var[:, 0], var[-1, 1:-1],
                              var[::-1, -1], var[0, ::-1][1:]])
        return func(lon), func(lat)

    x_rho = outputs.lon_rho

    ox_psi, oy_psi = outline(
        outputs.lon_psi[1:-1, 1:-1], outputs.lat_psi[1:-1, 1:-1])

    if min(x_rho.ravel()) < 180 and max(x_rho.ravel()) > 180:
        mid = 180
        x_rho = x_rho - 180
        ox_psi = ox_psi - 180
    else:
        mid = 0

    ax.plot(ox_psi, oy_psi, 'r', zorder=5)
    return ax

def plot_topo_bis(outputs, figure, ax):
    x_rho, y_rho = outputs.lon_rho, outputs.lat_rho

    axpos = ax.get_position()
    pos_x = axpos.x0
    pos_y = axpos.y0 - 0.08
    cax_width = axpos.width
    cax_height = 0.04

    ax.pcolormesh(x_rho, y_rho, np.ma.masked_where(outputs.mask_rho > 0, outputs.mask_rho),
                  vmin=0, vmax=1, zorder=2, cmap=plt.cm.copper_r)

    plotTopo = ax.pcolormesh(x_rho, y_rho, np.ma.masked_where(
        outputs.mask_rho < 1, outputs.h), zorder=2)

    pos_cax = figure.add_axes([pos_x, pos_y, cax_width, cax_height])
    cb = figure.colorbar(plotTopo, cax=pos_cax, orientation='horizontal')
    ax.grid(True, which='both', linewidth=2, color='gray', alpha=0.5, linestyle='--')

    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.1f}'))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.1f}'))

    return ax, cb

def plot_mask(outputs, figure, ax):
    x_rho, y_rho = outputs.lon_rho, outputs.lat_rho

    axpos = ax.get_position()
    pos_x = axpos.x0
    pos_y = axpos.y0 - 0.08
    cax_width = axpos.width
    cax_height = 0.04

    # Plotting the entire mask_rho without any condition
    plotMask = ax.pcolormesh(x_rho, y_rho, outputs.mask_rho, vmin=0, vmax=1, zorder=2, cmap=plt.cm.gray_r)

    pos_cax = figure.add_axes([pos_x, pos_y, cax_width, cax_height])
    cb = figure.colorbar(plotMask, cax=pos_cax, orientation='horizontal')
    
    ax.grid(True, which='both', linewidth=2, color='gray', alpha=0.5, linestyle='--')

    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.1f}'))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.1f}'))

    return ax, cb

#==================================================================================================

class ValueUpdater:
    def __init__(self, grid):
        # Flattening and stacking the longitude and latitude data
        x_rho, y_rho = grid.lon_rho, grid.lat_rho
        self.data = np.vstack((x_rho.flatten(), y_rho.flatten())).T
        self.x_rho = x_rho
        self.grid = grid
        self.plot_topo = plot_topo_bis
        self.plot_outline = plot_outline_ax
        
        plt.close('all')  # Close any existing figures

        # Create a new figure
        self.figure, self.ax = plt.subplots(figsize=(8, 8))
        
        # Initial plot
        self.ax = self.plot_outline(self.grid, self.ax)
        self.ax, self.cb = self.plot_topo(self.grid, self.figure, self.ax)
        self.ax.set_title("Bathymetry Edit")
        self.figure.canvas.mpl_connect('button_press_event', self.on_click)

        # Create and display widgets
        self.create_widgets()
        plt.show()
        
    def create_widgets(self):
        """Create and display widgets."""
        self.oldValue_widget = FloatText(value=0.0, description='Old value:', disabled=True)
        self.new_value_widget = FloatText(value=0.0, description='New value:')
        self.update_button = Button(description='Update Value')
        self.update_button.on_click(self.on_update_click)
        
        widget_box = VBox([self.oldValue_widget, self.new_value_widget, self.update_button])
        display(widget_box)
        
    def closest_node(self, node, nodes):
        """Find the closest node to the given point."""
        closest_index = distance.cdist([node], nodes).argmin()
        return closest_index, nodes[closest_index]
    
    def on_click(self, event):
        """Event handler for mouse clicks on the plot."""
        if event.inaxes == self.ax:
            pt = (event.xdata, event.ydata)
            prox_pt_id, prox_pt_value = self.closest_node(pt, self.data)
            p, mod = divmod(prox_pt_id, self.x_rho.shape[1])
            self.focusedPointX = p
            self.focusedPointY = mod
            self.oldValue_widget.value = self.grid.h[self.focusedPointX][self.focusedPointY]
    
    def on_update_click(self, b):
        """Handle the update button click event."""
        if hasattr(self, 'focusedPointX') and hasattr(self, 'focusedPointY'):
            new_value = self.new_value_widget.value
            self.grid.h[self.focusedPointX][self.focusedPointY] = new_value
            self.ax.clear()
            self.ax = self.plot_outline(self.grid, self.ax)
            self.ax, self.cb = self.plot_topo(self.grid, self.figure, self.ax)
            self.ax.set_title("Bathymetry Edit")
            self.figure.canvas.draw()  # Force a redraw of the figure
            self.figure.canvas.flush_events()  # Process any pending events

class RectangleSelectorEdition:
    def __init__(self, grid):
        # Flattening and stacking the longitude and latitude data
        x_rho, y_rho = grid.lon_rho, grid.lat_rho
        self.data = np.vstack((x_rho.flatten(), y_rho.flatten())).T
        self.x_rho = x_rho
        self.y_rho = y_rho
        self.grid = grid
        
        plt.close('all')  # Close any existing figures

        # Create a new figure
        self.figure, self.ax = plt.subplots(figsize=(8, 8))
        
        # Initial plot
        self.ax, self.cb = plot_topo_bis(self.grid, self.figure, self.ax)
        self.ax = plot_outline_ax(self.grid, self.ax)
        self.ax.set_title("Bathymetry Edit")

        # Set limits to fit the data
        self.ax.set_xlim(self.x_rho.min(), self.x_rho.max())
        self.ax.set_ylim(self.y_rho.min(), self.y_rho.max())
        
        # Create and display widgets
        self.create_widgets()
        
        # Add Rectangle Selector
        self.rect_selector = RectangleSelector(
            self.ax, 
            self.on_select,
            drawtype='box', 
            useblit=True,
            button=[1],  # Left mouse button
            minspanx=5, 
            minspany=5,
            spancoords='pixels',
            interactive=True
        )
        
        plt.show()

    def create_widgets(self):
        """Create and display widgets."""
        self.new_value_widget = FloatText(value=0.0, description='New value:')
        self.update_button = Button(description='Update Values')
        self.update_button.on_click(self.on_update_click)
        
        widget_box = VBox([self.new_value_widget, self.update_button])
        display(widget_box)

    def on_select(self, eclick, erelease):
        """Handle rectangle selection."""
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata

        # Ensure x1 < x2 and y1 < y2
        x_min, x_max = min(x1, x2), max(x1, x2)
        y_min, y_max = min(y1, y2), max(y1, y2)

        # Store selection bounds
        self.selection_bounds = (x_min, x_max, y_min, y_max)

    def on_update_click(self, b):
        """Handle the update button click event."""
        if hasattr(self, 'selection_bounds'):
            x_min, x_max, y_min, y_max = self.selection_bounds

            # Find indices within the rectangle
            mask_x = (self.x_rho >= x_min) & (self.x_rho <= x_max)
            mask_y = (self.y_rho >= y_min) & (self.y_rho <= y_max)
            mask = mask_x & mask_y

            # Update the grid values within the selected region
            self.grid.h[mask] = self.new_value_widget.value

            # Redraw the figure
            self.ax.clear()
            self.ax, self.cb = plot_topo_bis(self.grid, self.figure, self.ax)
            self.ax = plot_outline_ax(self.grid, self.ax)
            self.ax.set_title("Bathymetry Edit")
            self.figure.canvas.draw()  # Force a redraw of the figure

            # Clear selection bounds
            del self.selection_bounds

# -- Class for Mask edition --

class MaskUpdater:
    def __init__(self, grid):
        self.grid = grid
        x_rho, y_rho = grid.lon_rho, grid.lat_rho
        self.data = np.vstack((x_rho.flatten(), y_rho.flatten())).T
        self.x_rho = x_rho
        
        self.figure, self.ax = plt.subplots(figsize=(8, 8))
        self.plot_grid()

        self.oldValue_widget = FloatText(value=0.0, description='Old value:', disabled=True)
        self.new_value_widget = IntSlider(value=0, min=0, max=1, step=1, description='New value:')
        self.update_button = Button(description='Update Value')

        self.update_button.on_click(self.on_update_click)

        # Display widgets
        self.widget_box = VBox([self.oldValue_widget, self.new_value_widget, self.update_button])
        display(self.widget_box)
        
        self.figure.canvas.mpl_connect('button_press_event', self.on_click)

    def plot_grid(self):
        self.ax.clear()
        self.ax = plot_outline_ax(self.grid, self.ax)  # Assuming this function is defined
        self.ax, self.cb = plot_mask(self.grid, self.figure, self.ax)  # Assuming this function is defined
        self.ax.set_title("Mask Edit")
        self.figure.canvas.draw_idle()

    def closest_node(self, node, nodes):
        closest_index = distance.cdist([node], nodes).argmin()
        return closest_index, nodes[closest_index]

    def on_click(self, event):
        if event.inaxes == self.ax:
            pt = (event.xdata, event.ydata)
            prox_pt_id, prox_pt_value = self.closest_node(pt, self.data)
            p, mod = divmod(prox_pt_id, self.x_rho.shape[1])
            self.focusedPointX = p
            self.focusedPointY = mod
            self.oldValue_widget.value = self.grid.mask_rho[self.focusedPointX][self.focusedPointY]

    def on_update_click(self, b):
        if hasattr(self, 'focusedPointX') and hasattr(self, 'focusedPointY'):
            new_value = self.new_value_widget.value
            self.grid.mask_rho[self.focusedPointX][self.focusedPointY] = new_value
            self.plot_grid()
            self.figure.canvas.flush_events()