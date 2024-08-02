#--- Dependencies ---------------------------------------------------------
import configparser
import ipywidgets as widgets
from ipywidgets import VBox, Accordion, HTML
from IPython.display import display
import configparser

#--- Functions ------------------------------------------------------------
def load_namelist(filename):
    """
    Load the namelist variables with grid inputs config_grid.ini
    
    """
    config = configparser.ConfigParser()
    config.read(filename)
    
    # Grid section
    grid = config['Grid']
    tra_lon = float(grid['tra_lon'])
    tra_lat = float(grid['tra_lat'])
    size_x = int(grid['size_x'])
    size_y = int(grid['size_y'])
    nx = int(grid['nx'])
    ny = int(grid['ny'])
    rot = float(grid['rot'])
    
    # Smoothing section
    depth = config['Smoothing_params']
    hmin = float(depth['hmin'])
    hmax = float(depth['hmax'])
    interp_rad = int(depth['interp_rad'])
    rfact = float(depth['rfact'])
    smooth_meth = depth['smooth_meth']
    
    # Files section
    files = config['Files']
    topofile = files['topofile']
    shp_file = files['shp_file']
    output_file = files['output_file']
    
    # Single Connect section
    sgl_connect = config.getboolean('Single_Connect', 'sgl_connect')
    sgl_index1 = int(config['Single_Connect']['sgl_index1'])
    sgl_index2 = int(config['Single_Connect']['sgl_index2'])
    
    # Return all values in a dictionary
    return {
        'tra_lon': tra_lon,
        'tra_lat': tra_lat,
        'size_x': size_x,
        'size_y': size_y,
        'nx': nx,
        'ny': ny,
        'rot': rot,
        'hmin': hmin,
        'hmax': hmax,
        'interp_rad': interp_rad,
        'rfact': rfact,
        'smooth_meth': smooth_meth,
        'topofile': topofile,
        'shp_file': shp_file,
        'output_file': output_file,
        'sgl_connect': sgl_connect,
        'sgl_index1': sgl_index1,
        'sgl_index2': sgl_index2,
    }

#--- Widgets'stuff ---------------------------------------------------------

def setup_widgets():
    """
    Widgets for notebook make_grid.ipynb user's changes section
    """
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

    # Dictionary to store widget values after save button is clicked
    saved_config = {}

    # Define the save function
    def save_config(button):
        nonlocal saved_config  # Allow modification of saved_config inside this function
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
        }
        
        # Optionally, write to a config file
        config = configparser.ConfigParser()
        
        config['Grid'] = {
            'tra_lon': str(saved_config['tra_lon']),
            'tra_lat': str(saved_config['tra_lat']),
            'size_x': str(saved_config['size_x']),
            'size_y': str(saved_config['size_y']),
            'nx': str(saved_config['nx']),
            'ny': str(saved_config['ny']),
            'rot': str(saved_config['rot']),
        }
        
        config['Depth'] = {
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
        
        with open('config_grid_notebook.ini', 'w') as configfile:
            config.write(configfile)
        
        print("Configuration saved to 'config_grid_notebook.ini'")

    # Define the save button
    save_button = widgets.Button(description="Save")
    save_button.on_click(save_config)  # Attach save_config function to save_button click event

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
    accordion = Accordion(children=[
        grid_center_box, 
        grid_size_box, 
        grid_rotation_box, 
        depth_params_box, 
        file_params_box, 
        single_connect_box
    ])
    
    accordion.set_title(0, 'Grid Center')
    accordion.set_title(1, 'Grid Size')
    accordion.set_title(2, 'Grid Rotation')
    accordion.set_title(3, 'Smoothing Parameters')
    accordion.set_title(4, 'Files')
    accordion.set_title(5, 'Single Connect')

    # Display the Accordion and save button
    display(VBox([accordion, save_button]))
    
    # Return the function to get saved values
    return lambda: saved_config

# Appeler la fonction pour afficher les widgets
get_saved_config = setup_widgets()

# Utiliser la fonction get_saved_config pour obtenir les valeurs sauvegardées
saved_config = get_saved_config()
