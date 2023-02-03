# Widget
# cf. https://ipywidgets.readthedocs.io/en/stable/index.html
import ipywidgets as widgets
import make_grid_param as param
 
tra_lon = widgets.BoundedFloatText(
    value=param.tra_lon,
    min=-180,
    max=360,
    description='Longitude of the grid center:',
    disabled=False
)

tra_lat = widgets.BoundedFloatText(
    value=param.tra_lat,
    min=-90,
    max=90,
    description='Latitude of the grid center:',
    disabled=False
)

size_x = widgets.FloatText(
    value=param.size_x,
    step=1,
    description='X grid resolution [km]:',
    disabled=False
)

size_y = widgets.FloatText(
    value=param.size_y,
    step=1,
    description='Y grid resolution [km]:',
    disabled=False
)

nx = widgets.IntText(
    value=param.nx,
    step=1,
    description='nx:',
    disabled=False
)

ny = widgets.IntText(
    value=param.ny,
    step=1,
    description='ny:',
    disabled=False
)

rot = widgets.FloatText(
    value=param.rot,
    step=0.1,
    description='Rotation [degree]:',
    disabled=False
)


hmin = widgets.FloatText(
    value=param.hmin,
    step=0.1,
    description='Minimum depth [m]:',
    disabled=False
)


hmax = widgets.FloatText(
    value=param.hmax,
    step=0.1,
    description='Maximum depth [m]:',
    disabled=False
)

smth_rad = widgets.IntText(
    value=param.smth_rad,
    step=0.1,
    description='Smoothing radius [nb points] (usually between 2 and 8):',
    disabled=False
)

rfact = widgets.FloatText(
    value=param.rfact,
    step=0.2,
    description='Maximum r-fact to reach (the lower it is, the smoother it will be):',
    disabled=False
)

smooth_meth = widgets.Dropdown(
    options=['smooth', 'lsmooth', 'lsmooth_legacy', 'lsmooth2', 'lsmooth1', 'cond_rx0_topo'],
    value=param.smooth_meth,
    description='Smoothing method:',
    disabled=False,
)

topofile = widgets.Text(
    value=param.topofile,
    description='Topo file path:',
    disabled=False,
)

shp_file = widgets.Text(
    value=param.shp_file,
    description='Coastline file (for the mask) file path:',
    disabled=False,
)
