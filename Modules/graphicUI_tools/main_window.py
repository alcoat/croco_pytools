import numpy as np
import matplotlib
import matplotlib.pyplot
# We want matplotlib to use a wxPython backend
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from cartopy import crs as ccrs, feature as cfeature
import cartopy.io.shapereader as shpreader
# from threading import Thread
from time import sleep
import wx

from traits.api import *
from traitsui.api import View, Item,HGroup, Group, HSplit, Handler, FileEditor,DirectoryEditor,CheckListEditor
from traitsui.menu import NoButtons
from traitsui.wx.editor import Editor
#from traitsui.wx.basic_editor_factory import BasicEditorFactory
from traitsui.basic_editor_factory import BasicEditorFactory

from Modules.graphicUI_tools.button_actions import *
from Modules.tools_make_grid import EasyGrid, GetTopo, GetMask, topo_prt

from Modules.croco_class import CROCO

from Modules.graphicUI_tools.outputs_ui import OutputsUi
import make_grid_param as param

from Modules.map_tools.map_tools import plot_grid, plot_topo
        
#################################
#### FOR FIGURE #################
class _MPLFigureEditor(Editor):
    """
    from http://github.enthought.com/traitsui/tutorials/traits_ui_scientific_app.html
    """
    scrollable  = True

    def init(self, parent):
        self.control = self._create_canvas(parent)
        self.set_tooltip()

    def update_editor(self):
        pass

    def _create_canvas(self, parent):
        """ Create the MPL canvas. """
        # The panel lets us add additional controls.
        panel = wx.Panel(parent, -1, style=wx.CLIP_CHILDREN)
        sizer = wx.BoxSizer(wx.VERTICAL)
        panel.SetSizer(sizer)
        # matplotlib commands to create a canvas
        mpl_control = FigureCanvas(panel, -1, self.value)
        sizer.Add(mpl_control, 1, wx.LEFT | wx.TOP | wx.GROW)
        toolbar = NavigationToolbar2Wx(mpl_control)
        sizer.Add(toolbar, 0, wx.EXPAND)
        self.value.canvas.SetMinSize((10, 10))
        return panel

class MPLFigureEditor(BasicEditorFactory):
    """
    from http://github.enthought.com/traitsui/tutorials/traits_ui_scientific_app.html
    """
    klass = _MPLFigureEditor
#############################

class Inputs(HasTraits):
    """
    Inputs object
    """
    zview = Enum('grid outline', 'grid points', 'topo', '1/pm', '1/pn', 'angle', 'mask',
        desc="the data to view",
        label="View", )

    tra_lon = CFloat(param.tra_lon,
        desc="a central longitude",
        label="longitude", )

    tra_lat = CFloat(param.tra_lat,
        desc="a central latitude",
        label="latitude", )

    size_x = CFloat(param.size_x,
        desc="the mean distance along xi",
        label="x size [km]", )

    size_y = CFloat(param.size_y,
        desc="the mean distance along eta",
        label="y size [km]", )

    rot = CFloat(param.rot,
        desc="rotation about longitude, latitude",
        label="rotation [deg]", )

    nx = CInt(param.nx,
        desc="the number of points along xi",
        label="nx", )

    ny = CInt(param.ny,
        desc="the number of points along eta",
        label="ny", )

class Inputs_smth(HasTraits):
    """
    Inputs object for smoothing
    """
    depthmin=CFloat(param.hmin,
            desc="minimum depth",
            label="Minimum depth [m]",)

    depthmax=CFloat(param.hmax,
            desc="maximum depth",
            label="Maximum depth [m]",)

    smthr=CFloat(param.smth_rad,
            desc="smoothing radius",
            label="Smth radius [nb points]",)

    rfact=CFloat(param.rfact,
            desc="maximum r-factor",
            label="r-factor",)

    smooth= Enum('smooth', 'lsmooth', 'lsmooth_legacy', 'lsmooth2', 'lsmooth1', 'cond_rx0_topo',
        desc="smoothing method",
        label="Smoothing method", )

class Inputs_smth_c2c(HasTraits):
    """
    Inputs object for smoothing
    """
    smthr=CFloat(param.smth_rad,
            desc="smoothing radius",
            label="Smth radius [nb points]",)

    rfact=CFloat(param.rfact,
            desc="maximum r-factor",
            label="r-factor",)

    smooth= Enum('smooth', 'lsmooth', 'lsmooth_legacy', 'lsmooth2', 'lsmooth1', 'cond_rx0_topo',
        desc="smoothing method",
        label="Smoothing method", )



class Inputs_zm(HasTraits):
    """
    Inputs object
    """
    tra_lon = CFloat(param.zoom_tra_lon,
        desc="a central longitude",
        label="longitude", )

    tra_lat = CFloat(param.zoom_tra_lat,
        desc="a central latitude",
        label="latitude", )

    size_x = CFloat(param.zoom_size_x,
        desc="the mean distance along xi",
        label="x size [km]", )

    size_y = CFloat(param.zoom_size_y,
        desc="the mean distance along eta",
        label="y size [km]", )

    rot = CFloat(param.zoom_rot,
        desc="rotation about longitude, latitude",
        label="rotation [deg]", )

    nx = CInt(param.zoom_nx,
        desc="the number of points along xi",
        label="nx", )

    ny = CInt(param.zoom_ny,
        desc="the number of points along eta",
        label="ny", )


class Inputs_c2c(HasTraits):
    """
    Inputs object
    """
    coef = CInt(param.c2c_coef,
        desc="Refinement coefficient",
        label="Refinement coef")

    imin = CInt(param.c2c_imin,
        desc="Parent imin",
        label="imin", )

    imax = CInt(param.c2c_imax,
        desc="Parent imax",
        label="imax", )

    jmin = CInt(param.c2c_jmin,
        desc="Parent jmin",
        label="jmin", )

    jmax = CInt(param.c2c_jmax,
        desc="Parent jmax",
        label="jmax", )
#############################

class MainWindowHandler(Handler):
    def close(self, info, is_OK):
        if ( info.object.panel.compute_grid_thread
            and info.object.panel.compute_grid_thread.is_alive() ):
            info.object.panel.compute_grid_thread.wants_abort = True
            while info.object.panel.compute_grid_thread.is_alive():
                sleep(0.1)
            wx.Yield()

        return True

class MainWindow(HasTraits):
    """
    The main window, here go the instructions to create and destroy the application.
    """

    ##### Main instances #####
    figure = Instance(Figure)

    inputs = Instance(Inputs,())
    inputs_smth = Instance(Inputs_smth, ())
    inputs_zm = Instance(Inputs_zm, ())
    inputs_c2c = Instance(Inputs_c2c, ())
    inputs_smth_c2c = Instance(Inputs_smth_c2c, ())

    outputs = Instance(OutputsUi, ())
    opt_file = File(value=param.output_file,
                     label='Output file',
                     desc='Output path')
    opt_file_zm = File(value=param.zoom_output_file,
                     label='Output file',
                     desc='Output path')
    opt_dir = File(value=param.output_dir,
                     label='Output dir',
                     desc='Output path')

    results_string = String()

    # Button #
    compute_grid = Button("Compute grid") # Compute grid location
    compute_grid_thread = Instance(ComputeGridThread)


    compute_smooth = Button("Compute smoothing") # Compute grid smoothing
    compute_smooth_thread = Instance(ComputeSmthThread)

    compute_zm = Button("Compute child grid") # Compute offlinz zoom
    compute_zm_thread = Instance(ComputeZmThread)

    compute_c2c = Button("Compute child grid") # Compute AGRIF zoom
    compute_c2c_thread = Instance(ComputeC2cThread)
       # different save grid depending on the case #
    save_grid = Button("Save grid")
    save_grid_zm = Button("Save grid") # Offline zoom
    save_grid_c2c = Button("Save grid") # Agrif save button
    save_grid_thread = Instance(SaveGridThread)
    save_grid_zm_thread = Instance(SaveGridThread)
    save_grid_c2c_thread = Instance(SaveGridThread)



    # Get file #
    croco_file = File(value=param.output_file,
                     label='Parent grid file',
                     desc='Parent path and filename')

    shp_file=File(value=param.shp_file,
                     label='Coastline file',
                     desc='costline path')
    
    topo_file = File(value=param.topofile,
                     label='Topography file',
                     desc='path and topography filename')

    ## Functions and others ##
    easy = Instance(EasyGrid, ())
    get_topo = Instance(GetTopo, ())
    get_mask = Instance(GetMask, ())
    save2netcdf = Instance(CROCO, ())

    single_connect = Bool(value=True,label='Single connect')
    sglc_i = CInt(param.sgl_connect[1],label='i0')
    sglc_j = CInt(param.sgl_connect[2],label='j0')

    merge= CInt(param.zoom_merge, label='Merging area (nb points)')
    checklist = List(editor=CheckListEditor(values=['South','West', 'East', 'North'], cols=4 ))

    ##########################
    def _figure_default(self):
        figure = Figure() 
        #figure.add_axes([0.05, 0.04, 0.9, 0.92])
        return figure

#    def _panel_default(self):
#        return ControlPanel(figure=self.figure)


    panel = Group(
                 ######### First tab ########
                    Group(
                          Item(name='inputs', style='custom', show_label=False, springy=True),
                          '_',
                          Item(name='shp_file',editor=DirectoryEditor(entries=1),style='simple'),
                          '_',
                          Item(name='results_string', style='custom', show_label=False, springy=True, height=200 ),
                          '_',
                          Item(name='topo_file',editor=FileEditor(filter=['*.nc']),
                               style='simple', show_label=True, springy=True),
                          '_',
                          Item(name='compute_grid', show_label=False)
                          ,label="Configure grid", dock='tab'),
                  ########## Second tab ######
                    Group(
                          Item('inputs_smth',style='custom', show_label=False, springy=True),
                          '_',
                          HGroup('single_connect','sglc_i', 'sglc_j'),
                          '_',
                          Item(name='results_string', style='custom', show_label=False, springy=True, height=200 ),
                          '_',
                          Item(name='compute_smooth', show_label=False),
                          '_',
                          Item(name='opt_file',show_label=True),
                          Item(name='save_grid', show_label=False )
                          ,label="Topo smoothing",dock='tab'),
                  ########## Third tab ######
                    Group(
                          Item(name='croco_file',editor=FileEditor(filter=['*.nc']),
                               style='simple', show_label=True, springy=True),
#                          '_',
                          Item(name='inputs_zm', style='custom', show_label=False, springy=True),
#                          '_',
                          Item(name='shp_file',editor=DirectoryEditor(entries=1),style='simple'),
#                          '_',
                          Item(name='checklist',label='Open boundaries', style='custom'),
                          Item(name='merge',label='Merging area (nb points)',style='simple', springy=True),
#                          '_',
                          HGroup('single_connect','sglc_i', 'sglc_j'),
#                          '_',
                          Item(name='topo_file',editor=FileEditor(filter=['*.nc']),
                               style='simple', show_label=True, springy=True),
                          Item('inputs_smth',style='custom', show_label=False, springy=True),
#                          '_',
                          Item(name='compute_zm', show_label=False),
#                          '_',
                          Item(name='opt_file_zm',show_label=True),
                          Item(name='save_grid_zm', show_label=False ),
                          label="Create offline zoom",dock='tab'),
                  ########## Fourth tab ########
                    Group(
                          Item(name='croco_file',editor=FileEditor(filter=['*.nc*']),
                               style='simple', show_label=True, springy=True),
                          '_',
                          Item(name='inputs_c2c', style='custom', show_label=False, springy=True),
                          '_',
                          Item(name='shp_file',editor=DirectoryEditor(entries=1),style='simple'),
                          '_',
                          Item(name='checklist',label='Open boundaries', style='custom', id="custom"),
                          Item(name='merge',label='Merging area (nb points)',style='simple', springy=True),
                          '_',
                          HGroup('single_connect','sglc_i', 'sglc_j'),
                          '_',
                          Item(name='topo_file',editor=FileEditor(filter=['*.nc*']),
                               style='simple', show_label=True, springy=True),
                          Item('inputs_smth_c2c',style='custom', show_label=False, springy=True),
                          '_',
                          Item(name='compute_c2c', show_label=False),
                          '_',
                          Item(name='opt_dir',show_label=True),
                          Item(name='save_grid_c2c', show_label=False ),
                          label="Create zoom AGRIF",dock='tab'),
                    layout='tabbed') 
                    


    view = View(HSplit(Item('figure', editor=MPLFigureEditor(),dock='vertical',show_label=False),panel),
                resizable=True,
                height=0.90, width=0.75,
                buttons=NoButtons,
                title="Easygrid")


    def _compute_grid_fired(self):
        self.compute_grid_thread = ComputeGridThread()
        self.compute_grid_thread.inputs = self.inputs
        self.compute_grid_thread.outputs = self.outputs
        self.compute_grid_thread.outputs_file = self.opt_file
        if self.compute_grid_thread.inputs.zview=='topo':
            self.compute_grid_thread.topo_file = self.topo_file
            self.compute_grid_thread.shp_file = self.shp_file
            self.compute_grid_thread.topo = self.get_topo.topo
        elif self.compute_grid_thread.inputs.zview=='mask':
            self.compute_grid_thread.shp_file = self.shp_file
            self.compute_grid_thread.mask = self.get_mask.mask
        self.compute_grid_thread.display = self.add_line
        self.compute_grid_thread.easy = self.easy.easygrid
        self.compute_grid_thread.grid_show = self.grid_show
        self.compute_grid_thread.start()
        self.compute_grid_thread.join()
        self.compute_grid_thread.grid_show()

    def _compute_smooth_fired(self):
        self.compute_smooth_thread = ComputeSmthThread()
        self.compute_smooth_thread.inputs = self.inputs
        self.compute_smooth_thread.inputs_smth = self.inputs_smth
        self.compute_smooth_thread.outputs = self.outputs
        self.compute_smooth_thread.mask = self.get_mask.mask
        self.compute_smooth_thread.topo_file = self.topo_file
        self.compute_smooth_thread.topo = self.get_topo.topo
        self.compute_smooth_thread.shp_file = self.shp_file
        self.compute_smooth_thread.easy = self.easy.easygrid
        self.compute_smooth_thread.grid_show = self.grid_show
        self.compute_smooth_thread.single_connect = [self.single_connect,self.sglc_i,self.sglc_j]
        self.compute_smooth_thread.display = self.add_line
        self.compute_smooth_thread.start()
        self.compute_smooth_thread.join()
        self.compute_smooth_thread.grid_show(smooth=True)

    def _compute_zm_fired(self):
        self.compute_zm_thread = ComputeZmThread()
        self.compute_zm_thread.inputs = self.inputs_zm
        self.compute_zm_thread.inputs_smth = self.inputs_smth
        self.compute_zm_thread.outputs = self.outputs
        self.compute_zm_thread.mask = self.get_mask.mask
        self.compute_zm_thread.topo_file = self.topo_file
        self.compute_zm_thread.topo = self.get_topo.topo
        self.compute_zm_thread.match_topo = self.get_topo.match_topo
        self.compute_zm_thread.openb  = [self.checklist,self.merge]
        self.compute_zm_thread.shp_file = self.shp_file
        self.compute_zm_thread.easy = self.easy.easygrid
        self.compute_zm_thread.croco_file = self.croco_file
        self.compute_zm_thread.single_connect = [self.single_connect,self.sglc_i,self.sglc_j]
        self.compute_zm_thread.grid_show_zm = self.grid_show_zm
        self.compute_zm_thread.display = self.add_line
        self.compute_zm_thread.topo_prt = topo_prt(self.croco_file)
        prt_grd = topo_prt(self.croco_file)
        self.compute_zm_thread.start()
        self.compute_zm_thread.join()
        self.compute_zm_thread.grid_show_zm(prt_grd)

    def _compute_c2c_fired(self):
        self.compute_c2c_thread = ComputeC2cThread()
        self.compute_c2c_thread.inputs = self.inputs_c2c
        self.compute_c2c_thread.inputs_smth = self.inputs_smth_c2c
        self.compute_c2c_thread.outputs = self.outputs
        self.compute_c2c_thread.mask = self.get_mask.mask
        self.compute_c2c_thread.topo_file = self.topo_file
        self.compute_c2c_thread.topo = self.get_topo.topo
        self.compute_c2c_thread.match_topo = self.get_topo.match_topo
        self.compute_c2c_thread.shp_file = self.shp_file
        self.compute_c2c_thread.openb  = [self.checklist,self.merge]
        self.compute_c2c_thread.single_connect = [self.single_connect,self.sglc_i,self.sglc_j]
        self.compute_c2c_thread.nest = self.easy.AGRIFgrid
        self.compute_c2c_thread.topo_prt = topo_prt(self.croco_file)
        prt_grd = topo_prt(self.croco_file)
        self.compute_c2c_thread.grid_show = self.grid_show_zm
        self.compute_c2c_thread.display = self.add_line
        self.compute_c2c_thread.start()
        self.compute_c2c_thread.join()
        self.compute_c2c_thread.grid_show(prt_grd)

    def _save_grid_fired(self):
        self.save_grid_thread = SaveGridThread()
        self.save_grid_thread.inputs = self.inputs
        self.save_grid_thread.outputs = self.outputs
        self.save_grid_thread.outputs_file = self.opt_file
        self.save_grid_thread.display = self.add_line
        self.save_grid_thread.prt_grd= None
        self.save_grid_thread.save2netcdf = self.save2netcdf
        self.save_grid_thread.start()

    def _save_grid_zm_fired(self):
        self.save_grid_zm_thread = SaveGridThread()
        self.save_grid_zm_thread.inputs = self.inputs_zm
        self.save_grid_zm_thread.outputs = self.outputs
        self.save_grid_zm_thread.outputs_file = self.opt_file_zm
        self.save_grid_zm_thread.display = self.add_line
        self.save_grid_zm_thread.prt_grd=[False]
        self.save_grid_zm_thread.save2netcdf = self.save2netcdf
        self.save_grid_zm_thread.start()

    def _save_grid_c2c_fired(self):
        self.save_grid_c2c_thread = SaveGridThread()
        self.save_grid_c2c_thread.inputs = self.inputs
        self.save_grid_c2c_thread.outputs = self.outputs
        self.save_grid_c2c_thread.outputs_file = self.opt_dir
        self.save_grid_c2c_thread.display = self.add_line
        self.save_grid_c2c_thread.prt_grd=[True,self.croco_file,self.inputs_c2c.coef,self.inputs_c2c.imin
                ,self.inputs_c2c.imax,self.inputs_c2c.jmin,self.inputs_c2c.jmax]
        self.save_grid_c2c_thread.save2netcdf = self.save2netcdf
        self.save_grid_c2c_thread.start()


######################################
######## END Saving Thread ###########
######################################

######################################
######## START display Threads #######
######################################

    def add_line(self, string):
        """
        Adds a line to the textbox display.
        """
        self.results_string = (string + "\n" + self.results_string)[0:1000]

    def grid_show(self,smooth=None):
        self.figure.clf()
        if smooth is None:
            plot_grid(self.outputs, self.figure, self.inputs.zview, self.shp_file)
        else:
            plot_topo(self.outputs, self.figure)
        
        wx.CallAfter(self.figure.canvas.draw)

    def grid_show_zm(self,prt_grd):
        def outline(lon, lat):
            '''
            Return lon, lat of perimeter around the grid
            '''
            def func(var):
                return np.hstack([var[:, 0], var[-1, 1:-1],
                                  var[::-1, -1], var[0, ::-1][1:]])
            return func(lon), func(lat)

        prt_xr, prt_yr = prt_grd.lon_rho,prt_grd.lat_rho
        prt_xp, prt_yp = prt_grd.lon_psi,prt_grd.lat_psi

        chd_xr, chd_yr = self.outputs.lon_rho,self.outputs.lat_rho
        chd_xp, chd_yp = self.outputs.lon_psi,self.outputs.lat_psi

        prt_oxp, prt_oyp  = outline(prt_grd.lon_psi[1:-1, 1:-1], prt_grd.lat_psi[1:-1, 1:-1])
        chd_oxp, chd_oyp  = outline(self.outputs.lon_psi[1:-1, 1:-1], self.outputs.lat_psi[1:-1, 1:-1])

        if min(prt_xr.ravel())<180 and max(prt_xr.ravel())>180:
            # check if grid cross lon=180 and change central lon for cartopy
            mid=180
            prt_xr=prt_xr-180;chd_xr=chd_xr-180
            prt_oxp=prt_oxp-180;chd_oxp=chd_oxp-180
        else:
            mid=0


        self.figure.clf()

        ax=self.figure.add_axes(rect=[0.1, 0.04, 0.80, 0.9],projection=ccrs.PlateCarree())

        ax.plot(prt_oxp, prt_oyp, 'r',zorder=5) #, ax=self.figure.gca(), zorder=3)
        ax.plot(chd_oxp, chd_oyp, 'r',zorder=5)
        # Colorbar position
        axpos = ax.get_position()
        pos_x = axpos.x0 # + 0.25*axpos.width
        pos_y = axpos.y0-0.08 #-axpos.height - 0.02
        cax_width = axpos.width
        cax_height = 0.04

        ax.pcolormesh(prt_xr, prt_yr,np.ma.masked_where(prt_grd.mask_rho>0, prt_grd.mask_rho) \
                ,vmin=0,vmax=1, zorder=2, cmap=matplotlib.pyplot.cm.copper_r)
        ax.pcolormesh(prt_xr, prt_yr, np.ma.masked_where(prt_grd.mask_rho<1, prt_grd.h) \
                ,vmin=np.nanmin(prt_grd.h),vmax=np.nanmax(prt_grd.h), zorder=3 )

        ax.pcolormesh(chd_xr, chd_yr,np.ma.masked_where(self.outputs.mask_rho>0, self.outputs.mask_rho) \
                ,vmin=0,vmax=1, zorder=3, cmap=matplotlib.pyplot.cm.copper_r)

        cb = ax.pcolormesh(chd_xr, chd_yr, np.ma.masked_where(self.outputs.mask_rho<1, self.outputs.h) \
                ,vmin=np.nanmin(prt_grd.h),vmax=np.nanmax(prt_grd.h),zorder=4 )

        pos_cax = self.figure.add_axes([pos_x,pos_y,cax_width,cax_height])
        self.figure.colorbar(cb,cax=pos_cax,orientation='horizontal')
        gl=ax.gridlines(draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False

        wx.CallAfter(self.figure.canvas.draw)

######################

