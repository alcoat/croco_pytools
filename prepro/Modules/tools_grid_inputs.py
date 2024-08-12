#--- Dependencies ---------------------------------------------------------
import configparser
import ipywidgets as widgets
from ipywidgets import VBox, Accordion, HTML
from IPython.display import display
import configparser

#--- Functions ------------------------------------------------------------
import os

def load_namelist(filename):
    """
    Load the namelist variables from the specified configuration file.
    
    Parameters:
    - filename (str): Path to the configuration file.
    
    Returns:
    - dict: A dictionary with the configuration parameters.
    """
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"The file '{filename}' does not exist.")
    
    config = configparser.ConfigParser()
    config.read(filename)
    
    # Initialize the return dictionary
    params = {}
    
    # Load Grid section
    if 'Grid' in config:
        grid = config['Grid']
        try:
            params['tra_lon'] = float(grid.get('tra_lon', 0))
            params['tra_lat'] = float(grid.get('tra_lat', 0))
            params['size_x'] = int(grid.get('size_x', 0))
            params['size_y'] = int(grid.get('size_y', 0))
            params['nx'] = int(grid.get('nx', 0))
            params['ny'] = int(grid.get('ny', 0))
            params['rot'] = float(grid.get('rot', 0))
        except ValueError as e:
            raise ValueError(f"Error parsing 'Grid' section: {e}")
    
    # Load Smoothing section
    if 'Smoothing_params' in config:
        depth = config['Smoothing_params']
        try:
            params['hmin'] = float(depth.get('hmin', 0))
            params['hmax'] = float(depth.get('hmax', 0))
            params['interp_rad'] = int(depth.get('interp_rad', 0))
            params['rfact'] = float(depth.get('rfact', 0))
            params['smooth_meth'] = depth.get('smooth_meth', 'default_method')
        except ValueError as e:
            raise ValueError(f"Error parsing 'Smoothing_params' section: {e}")
    
    # Load Files section
    if 'Files' in config:
        files = config['Files']
        params['topofile'] = files.get('topofile', '')
        params['shp_file'] = files.get('shp_file', '')
        params['output_file'] = files.get('output_file', '')
    
    # Load Single Connect section
    if 'Single_Connect' in config:
        single_connect = config['Single_Connect']
        try:
            params['sgl_connect'] = single_connect.getboolean('sgl_connect')
            params['sgl_index1'] = int(single_connect.get('sgl_index1', 0))
            params['sgl_index2'] = int(single_connect.get('sgl_index2', 0))
        except ValueError as e:
            raise ValueError(f"Error parsing 'Single_Connect' section: {e}")
    
    return params

#--- Widgets'stuff ---------------------------------------------------------

from ipywidgets import VBox, HBox, HTML, Accordion
from IPython.display import clear_output

from ipywidgets import widgets, VBox, HBox, HTML, Accordion, Button, Output
from IPython.display import display, clear_output
import configparser

def setup_widgets(parent_grid=False):
    """
    Widgets for notebook make_grid.ipynb user's changes section
    """
    # Initialize the saved_config dictionary
    saved_config = {}
    
    # Define the widgets for each section
    tra_lon = widgets.FloatText(description='Longitude:', value=15)
    tra_lat = widgets.FloatText(description='Latitude:', value=-32)
    size_x = widgets.IntText(description='Size X (km):', value=1556)
    size_y = widgets.IntText(description='Size Y (km):', value=1334)
    nx = widgets.IntText(description='NX:', value=39)
    ny = widgets.IntText(description='NY:', value=40)
    rot = widgets.FloatText(description='Rotation:', value=0)
    hmin = widgets.FloatText(description='Min Depth (m):', value=50)
    hmax = widgets.FloatText(description='Max Depth (m):', value=6000)
    interp_rad = widgets.IntText(description='Interp Radius:', value=2)
    rfact = widgets.FloatText(description='R-Fact:', value=0.2)
    smooth_meth = widgets.Dropdown(
        description='Smooth Method:',
        options=['smooth', 'lsmooth', 'lsmooth_legacy', 'lsmooth2', 'lsmooth1', 'cond_rx0_topo'],
        value='lsmooth'
    )
    topofile = widgets.Text(description='Topo File:', value='../../DATASETS_CROCOTOOLS/Topo/etopo2.nc')
    shp_file = widgets.Text(description='Shapefile:', value='../../DATASETS_CROCOTOOLS/gshhs/GSHHS_shp/i/GSHHS_i_L1.shp')
    output_file = widgets.Text(description='Output File:', value='../../CROCO_FILES/croco_grd.nc')
    sgl_connect = widgets.Checkbox(description='Single Connect?', value=False)
    sgl_index1 = widgets.IntText(description='SGI Index 1:', value=20)
    sgl_index2 = widgets.IntText(description='SGI Index 2:', value=20)

    # Create a button to save the settings
    save_button = widgets.Button(description="Save Settings")
    output = widgets.Output()

    # Function to handle saving values
    def save_values(b):
        nonlocal saved_config  # Allow modification of saved_config inside this function

        if parent_grid:
            # Update settings with current values
            settings = {
                'NORTH': north_checkbox.value,
                'SOUTH': south_checkbox.value,
                'WEST': west_checkbox.value,
                'EAST': east_checkbox.value,
                'MERGING_AREA': merging_area_text.value
            }
            # Convert settings to directions
            directions = []
            if settings['WEST']:
                directions.append('West')
            if settings['EAST']:
                directions.append('East')
            if settings['SOUTH']:
                directions.append('South')
            if settings['NORTH']:
                directions.append('North')
            settings['dirs'] = directions
        else:
            settings = {}

        # Save to saved_config
        saved_config = {
            'tra_lon': tra_lon.value,
            'tra_lat': tra_lat.value,
            'size_x': size_x.value,
            'size_y': size_y.value,
            'nx': nx.value,
            'ny': ny.value,
            'rot': rot.value,
            'hmin': hmin.value,
            'hmax': hmax.value,
            'interp_rad': interp_rad.value,
            'rfact': rfact.value,
            'smooth_meth': smooth_meth.value,
            'topofile': topofile.value,
            'shp_file': shp_file.value,
            'output_file': output_file.value,
            'sgl_connect': sgl_connect.value,
            'sgl_index1': sgl_index1.value,
            'sgl_index2': sgl_index2.value,
            'open_boundary': settings,
        }
        
        if parent_grid:
            saved_config['parent_grid'] = parent_grid_field.value
        
        # Display confirmation message
        with output:
            clear_output()
            print("Settings have been saved.")
            print(saved_config)  # Display current settings for verification

    save_button.on_click(save_values)

    # Simplified "Open Boundary Conditions" section if parent_grid=True
    open_boundary_box = None
    if parent_grid:
        parent_grid_field = widgets.Text(description='Parent Grid:', value='../../CROCO_FILES/croco_grd.nc')
        output_file = widgets.Text(description='Output File:', value='../../CROCO_FILES/croco_chd_grd.nc')

        # Settings for the open boundary conditions
        settings = {
            'NORTH': True,
            'SOUTH': True,
            'WEST': True,
            'EAST': True,
            'MERGING_AREA': 5  # Set default value to 5
        }
    
        # Create checkboxes for each border
        north_checkbox = widgets.Checkbox(value=settings['NORTH'], description='North Border ⬆️')
        south_checkbox = widgets.Checkbox(value=settings['SOUTH'], description='South Border ⬇️')
        west_checkbox = widgets.Checkbox(value=settings['WEST'], description='West Border ⬅️')
        east_checkbox = widgets.Checkbox(value=settings['EAST'], description='East Border ➡️')
    
        # Create a field for the "Merging Area" parameter as an integer with default value 5
        merging_area_text = widgets.IntText(
            value=settings['MERGING_AREA'], 
            layout=widgets.Layout(width='100px')  # Adjust width of the input field
        )
        merging_area_label = widgets.Label(
            value='Merging Area:', 
            layout=widgets.Layout(width='150px', display='flex', justify_content='flex-end')  # Adjust label width and alignment
        )
    
        # Create an HBox for the label and input field
        merging_area_box = widgets.HBox(children=[
            merging_area_label, 
            merging_area_text
        ])

        # Create open boundary conditions box
        open_boundary_box = VBox([
            north_checkbox, 
            south_checkbox, 
            west_checkbox, 
            east_checkbox,
            merging_area_box
        ])
        open_boundary_box.layout.border = ''
        open_boundary_box.layout.padding = ''
        open_boundary_box.layout.border_radius = ''
        open_boundary_box.layout.margin = '0px'
    
    # Define sections with HTML for separation and titles
    def create_section(title, *widgets):
        section = VBox([HTML(f"<h3>{title}</h3>")] + list(widgets))
        return section

    # Create sections
    grid_center_box = create_section("Grid Center", tra_lon, tra_lat)
    grid_size_box = create_section("Grid Size", size_x, size_y, nx, ny)
    grid_rotation_box = create_section("Grid Rotation", rot)
    depth_params_box = create_section("Smoothing Parameters", hmin, hmax, interp_rad, rfact, smooth_meth)
    
    file_params_box = create_section("Files", topofile, shp_file, output_file)
    
    single_connect_box = create_section("Single Connect", sgl_connect, sgl_index1, sgl_index2)

    # Create the Accordion widget
    children = [grid_center_box, grid_size_box, grid_rotation_box, depth_params_box, file_params_box]
    if open_boundary_box:
        children.append(create_section("Open Boundary Conditions", open_boundary_box))
    children.append(single_connect_box)

    accordion = Accordion(children=children)
    
    accordion.set_title(0, 'Grid Center')
    accordion.set_title(1, 'Grid Size')
    accordion.set_title(2, 'Grid Rotation')
    accordion.set_title(3, 'Smoothing Parameters')
    accordion.set_title(4, 'Files')
    if open_boundary_box:
        accordion.set_title(len(children) - 2, 'Open Boundary Conditions')
    accordion.set_title(len(children) - 1, 'Single Connect')

    # Display the Accordion and save button
    display(VBox([accordion, save_button, output]))
    
    # Return the function to get saved values
    return lambda: saved_config

# - - - - - - - - - - - - - - - - - - - - - - -

def setup_widgets_agrif(parent_grid=False):
    """
    Widgets for notebook make_grid.ipynb user's changes section with AGRIF parameters.
    """
    # Initialize the saved_config dictionary
    saved_config = {}
    
    # Define the AGRIF parameters widgets
    coef = widgets.IntText(description='Refinement Coefficient:', value=3)
    imin = widgets.IntText(description='Imin:', value=10)
    imax = widgets.IntText(description='Imax:', value=27)
    jmin = widgets.IntText(description='Jmin:', value=8)
    jmax = widgets.IntText(description='Jmax:', value=28)

    # Define other widgets that are still relevant
    hmin = widgets.FloatText(description='Min Depth (m):', value=50)
    hmax = widgets.FloatText(description='Max Depth (m):', value=6000)
    interp_rad = widgets.IntText(description='Interp Radius:', value=2)
    rfact = widgets.FloatText(description='R-Fact:', value=0.2)
    smooth_meth = widgets.Dropdown(
        description='Smooth Method:',
        options=['smooth', 'lsmooth', 'lsmooth_legacy', 'lsmooth2', 'lsmooth1', 'cond_rx0_topo'],
        value='lsmooth'
    )
    topofile = widgets.Text(description='Topo File:', value='../../DATASETS_CROCOTOOLS/Topo/etopo2.nc')
    shp_file = widgets.Text(description='Shapefile:', value='../../DATASETS_CROCOTOOLS/gshhs/GSHHS_shp/i/GSHHS_i_L1.shp')
    output_file = widgets.Text(description='Output File:', value='../../CROCO_FILES/croco_grd.nc')
    sgl_connect = widgets.Checkbox(description='Single Connect?', value=False)
    sgl_index1 = widgets.IntText(description='SGI Index 1:', value=20)
    sgl_index2 = widgets.IntText(description='SGI Index 2:', value=20)

    # Optional parent grid text field
    if parent_grid:
        parent_grid_field = widgets.Text(description='Parent Grid:', value='../../CROCO_FILES/croco_grd.nc')
        output_file = widgets.Text(description='Output File:', value='../../CROCO_FILES/croco_chd_grd.nc')

    # Settings for the open boundary conditions
    settings = {
        'NORTH': True,
        'SOUTH': True,
        'WEST': True,
        'EAST': True,
        'MERGING_AREA': 5  # Set default value to 5
    }

    # Create checkboxes for each border
    north_checkbox = widgets.Checkbox(value=settings['NORTH'], description='North Border ⬆️')
    south_checkbox = widgets.Checkbox(value=settings['SOUTH'], description='South Border ⬇️')
    west_checkbox = widgets.Checkbox(value=settings['WEST'], description='West Border ⬅️')
    east_checkbox = widgets.Checkbox(value=settings['EAST'], description='East Border ➡️')

    # Create a field for the "Merging Area" parameter as an integer with default value 5
    merging_area_text = widgets.IntText(
        value=settings['MERGING_AREA'], 
        layout=widgets.Layout(width='100px')  # Adjust width of the input field
    )
    merging_area_label = widgets.Label(
        value='Merging Area:', 
        layout=widgets.Layout(width='150px', display='flex', justify_content='flex-end')  # Adjust label width and alignment
    )

    # Create an HBox for the label and input field
    merging_area_box = widgets.HBox(children=[
        merging_area_label, 
        merging_area_text
    ])

    # Create a button to save the settings
    save_button = widgets.Button(description="Save Settings")

    # Function to handle saving values
    def save_values(b):
        nonlocal saved_config  # Allow modification of saved_config inside this function

        # Update settings with current values
        settings['NORTH'] = north_checkbox.value
        settings['SOUTH'] = south_checkbox.value
        settings['WEST'] = west_checkbox.value
        settings['EAST'] = east_checkbox.value
        settings['MERGING_AREA'] = merging_area_text.value
        
        # Convert settings to directions
        directions = []
        if settings['WEST']:
            directions.append('West')
        if settings['EAST']:
            directions.append('East')
        if settings['SOUTH']:
            directions.append('South')
        if settings['NORTH']:
            directions.append('North')
        settings['dirs'] = directions
        
        # Save to saved_config
        saved_config = {
            'coef': coef.value,
            'imin': imin.value,
            'imax': imax.value,
            'jmin': jmin.value,
            'jmax': jmax.value,
            'hmin': hmin.value,
            'hmax': hmax.value,
            'interp_rad': interp_rad.value,
            'rfact': rfact.value,
            'smooth_meth': smooth_meth.value,
            'topofile': topofile.value,
            'shp_file': shp_file.value,
            'output_file': output_file.value,
            'sgl_connect': sgl_connect.value,
            'sgl_index1': sgl_index1.value,
            'sgl_index2': sgl_index2.value,
            'open_boundary': settings,
        }
        
        if parent_grid:
            saved_config['parent_grid'] = parent_grid_field.value
        
        # Display confirmation message
        with output:
            clear_output()
            print("Settings have been saved.")
            print(saved_config)  # Display current settings for verification

    save_button.on_click(save_values)

    # Simplified "Open Boundary Conditions" section if parent_grid=True
    if parent_grid:
        open_boundary_box = VBox([
            north_checkbox, 
            south_checkbox, 
            west_checkbox, 
            east_checkbox,
            merging_area_box
        ])
        open_boundary_box.layout.border = ''
        open_boundary_box.layout.padding = ''
        open_boundary_box.layout.border_radius = ''
        open_boundary_box.layout.margin = '0px'
        # Do not add save button here; use the global save button
    else:
        # Create a box around the widgets with custom styling
        border_box = widgets.VBox(children=[
            north_checkbox, 
            south_checkbox, 
            west_checkbox, 
            east_checkbox,
            merging_area_box,
            save_button
        ])
        
        # Apply CSS styles to the box
        border_box.layout.border = '2px solid #007bff'
        border_box.layout.padding = '10px'
        border_box.layout.border_radius = '5px'
        border_box.layout.margin = '10px'
        border_box.layout.width = '400px'  # Increase width to ensure full visibility

        # Adjust width of individual widgets
        north_checkbox.layout.width = '100%'
        south_checkbox.layout.width = '100%'
        west_checkbox.layout.width = '100%'
        east_checkbox.layout.width = '100%'
        merging_area_text.layout.width = '100px'  # Set a specific width for the input field

        # Create an output area for displaying messages
        output = widgets.Output()

        # Add save button to the border_box
        border_box.children += (save_button,)

    # Define the save function
    def save_config(button):
        nonlocal saved_config  # Allow modification of saved_config inside this function
        saved_config = {
            'coef': coef.value,
            'imin': imin.value,
            'imax': imax.value,
            'jmin': jmin.value,
            'jmax': jmax.value,
            'hmin': hmin.value,
            'hmax': hmax.value,
            'interp_rad': interp_rad.value,
            'rfact': rfact.value,
            'smooth_meth': smooth_meth.value,
            'topofile': topofile.value,
            'shp_file': shp_file.value,
            'output_file': output_file.value,
            'sgl_connect': sgl_connect.value,
            'sgl_index1': sgl_index1.value,
            'sgl_index2': sgl_index2.value,
        }
        
        if parent_grid:
            saved_config['parent_grid'] = parent_grid_field.value
        
        # Optionally, write to a config file
        config = configparser.ConfigParser()
        
        config['AGRIF'] = {
            'coef': str(saved_config['coef']),
            'imin': str(saved_config['imin']),
            'imax': str(saved_config['imax']),
            'jmin': str(saved_config['jmin']),
            'jmax': str(saved_config['jmax']),
        }
        
        config['Smoothing_params'] = {
            'hmin': str(saved_config['hmin']),
            'hmax': str(saved_config['hmax']),
            'interp_rad': str(saved_config['interp_rad']),
            'rfact': str(saved_config['rfact']),
            'smooth_meth': saved_config['smooth_meth'],
        }
        
        config['Files'] = {
            'topofile': saved_config['topofile'],
            'shp_file': saved_config['shp_file'],
            'output_file': saved_config['output_file'],
        }
        
        config['Single_Connect'] = {
            'sgl_connect': str(saved_config['sgl_connect']),
            'sgl_index1': str(saved_config['sgl_index1']),
            'sgl_index2': str(saved_config['sgl_index2']),
        }

        if parent_grid:
            config['Parent_Grid'] = {
                'parent_grid': saved_config['parent_grid'],
            }
        
        with open('config_grid_agrif_notebook.ini', 'w') as configfile:
            config.write(configfile)
        
        print("Configuration saved to 'config_grid_notebook.ini'")

    # Define sections with HTML for separation and titles
    def create_section(title, *widgets):
        section = VBox([HTML(f"<h3>{title}</h3>")] + list(widgets))
        return section

    # Create sections
    agrif_params_box = create_section("AGRIF Parameters", coef, imin, imax, jmin, jmax)
    depth_params_box = create_section("Smoothing Parameters", hmin, hmax, interp_rad, rfact, smooth_meth)
    
    if parent_grid:
        file_params_box = create_section("Files", topofile, shp_file, output_file, parent_grid_field)
        open_boundary_box = create_section("Open Boundary Conditions", open_boundary_box)
    else:
        file_params_box = create_section("Files", topofile, shp_file, output_file)
        open_boundary_box = create_section("Open Boundary Conditions", border_box)
    
    single_connect_box = create_section("Single Connect", sgl_connect, sgl_index1, sgl_index2)

    # Create the Accordion widget
    children = [agrif_params_box, depth_params_box, file_params_box]
    if open_boundary_box:
        children.append(open_boundary_box)  # Append Open Boundary Conditions if present
    children.append(single_connect_box)

    accordion = Accordion(children=children)
    
    accordion.set_title(0, 'AGRIF Parameters')
    accordion.set_title(1, 'Smoothing Parameters')
    accordion.set_title(2, 'Files')
    if open_boundary_box:
        accordion.set_title(3, 'Open Boundary Conditions')
    accordion.set_title(len(children) - 1, 'Single Connect')

    # Display the Accordion and save button
    display(VBox([accordion, save_button]))
    
    # Return the function to get saved values
    return lambda: saved_config
