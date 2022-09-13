import numpy as np
import matplotlib
# We want matplotlib to use a wxPython backend
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_wx import NavigationToolbar2Wx

from threading import Thread
from time import sleep
import wx

from traits.api import *
from traitsui.api import View, Item, Group, HSplit, Handler, EnumEditor, FileEditor,DirectoryEditor
from traitsui.menu import NoButtons
from traitsui.wx.editor import Editor
from traitsui.basic_editor_factory import BasicEditorFactory

import matplotlib.pyplot as plt
import tools
import tools_topo

class ComputeGridThread(Thread):
    """
    This is the worker thread that 
    spawns the easygrid processing job
    """
    wants_abort = False

    def run(self):
        """
        Runs the easygrid computing loop
        """
        self.bmapoptions
        self.easy(self.inputs, self.outputs)
        if self.inputs.zview == 'topo':
            self.topo(self.outputs, self.topo_file)
        if self.inputs.zview == 'mask':
            self.mask(self.outputs,self.bmapoptions,self.gshhs_file)

        dx = (1 / self.outputs.pm.max(), 1 / self.outputs.pm.min())
        dy = (1 / self.outputs.pn.max(), 1 / self.outputs.pn.min())
        easyparam = (self.inputs.nx,      self.inputs.ny,
                     self.inputs.tra_lon, self.inputs.tra_lat,
                     self.inputs.size_x,  self.inputs.size_y,
                     self.inputs.rot)

        self.display('----- min=dy %.2f,  max=dy %.2f\n' % (dy))
        self.display('----- min=dx %.2f,  max=dx %.2f'   % (dx))
        self.display(''.join(('--- computing grid (nx=%i, ny=%i, ',
                              'lon=%.2f, lat=%.2f, sx=%.2f, sy=%.2f, ',
                              'rot=%.2f)')) % easyparam)

class ComputeSmthThread(Thread):
    """
    This is the worker thread that
    spawns the smoothing processing job
    """
    wants_abort = False

    def run(self):
        """
        Runs the smoothing computing loop
        """
        self.bmapoptions
        self.easy(self.inputs, self.outputs)
        self.mask(self.outputs,self.bmapoptions,self.gshhs_file,sgl_connect=self.single_connect)
        self.topo(self.outputs, self.topo_file,smooth=self.inputs_smth)

        dx = (1 / self.outputs.pm.max(), 1 / self.outputs.pm.min())
        dy = (1 / self.outputs.pn.max(), 1 / self.outputs.pn.min())

        easyparam = (self.inputs.nx,      self.inputs.ny,
                     self.inputs.tra_lon, self.inputs.tra_lat,
                     self.inputs.size_x,  self.inputs.size_y,
                     self.inputs.rot)

        self.display('----- min=dy %.2f,  max=dy %.2f\n' % (dy))
        self.display('----- min=dx %.2f,  max=dx %.2f'   % (dx))
        self.display(''.join(('--- computing grid (nx=%i, ny=%i, ',
                              'lon=%.2f, lat=%.2f, sx=%.2f, sy=%.2f, ',
                              'rot=%.2f)')) % easyparam)

class ComputeZmThread(Thread):
    """
    This is the worker thread that 
    spawns the easygrid processing job
    """
    wants_abort = False

    def run(self):
        """
        Runs the easygrid computing loop
        """
 
#        self.prt_grd=tools_topo.topo_prt(self.croco_file) 
        self.bmapoptions
        self.easy(self.inputs, self.outputs)
        self.mask(self.outputs,self.bmapoptions,self.gshhs_file,sgl_connect=self.single_connect)
        self.topo(self.outputs, self.topo_file,smooth=self.inputs_smth)

        dx = (1 / self.outputs.pm.max(), 1 / self.outputs.pm.min())
        dy = (1 / self.outputs.pn.max(), 1 / self.outputs.pn.min())
        easyparam = (self.inputs.nx,      self.inputs.ny,
                     self.inputs.tra_lon, self.inputs.tra_lat,
                     self.inputs.size_x,  self.inputs.size_y,
                     self.inputs.rot)

        self.display('----- min=dy %.2f,  max=dy %.2f\n' % (dy))
        self.display('----- min=dx %.2f,  max=dx %.2f'   % (dx))
        self.display(''.join(('--- computing grid (nx=%i, ny=%i, ',
                              'lon=%.2f, lat=%.2f, sx=%.2f, sy=%.2f, ',
                              'rot=%.2f)')) % easyparam)

class ComputeC2cThread(Thread):
    """
    This is the worker thread that
    spawns the smoothing processing job
    """
    wants_abort = False

    def run(self):
        """
        Runs the smoothing computing loop
        """
#        if (self.inputs.nx-1)%self.inputs.coef != 0 or (self.inputs.ny-1)%self.inputs.coef != 0 :
#            if (self.inputs.nx-1)%self.inputs.coef!=0 :
#                self.display('--- Error: nx need to be set as follow: \n')
#                self.display(' nx = N*refine_coef + 1 --with N an integer \n')
#            if (self.inputs.ny-1)%self.inputs.coef != 0 :
#                self.display('--- Error: nx need to be set as follow: \n')
#                self.display(' ny = N*refine_coef + 1 -- with N an integer \n')
#        else:

        self.bmapoptions #read coastline option
        self.nest(self.topo_prt,self.inputs,self.outputs)
        self.mask(self.outputs,self.bmapoptions,self.gshhs_file,sgl_connect=self.single_connect)
        self.topo(self.outputs, self.topo_file,smooth=self.inputs_smth,hmin=np.nanmin(self.topo_prt.h),hmax=np.nanmax(self.topo_prt.h))
        self.match_topo(self.topo_prt,self.outputs,self.openb) 

        dx = (1 / self.outputs.pm.max(), 1 / self.outputs.pm.min())
        dy = (1 / self.outputs.pn.max(), 1 / self.outputs.pn.min())


class SaveGridThread(Thread):
    """
    
    """
    wants_abort = False

    def run(self):
        """
        Runs the easygrid save to nc loop
        """
        self.save2netcdf(self.outputs_dir,self.inputs, self.outputs,prt_grd=self.prt_grd)

        if self.prt_grd is None or self.prt_grd[0]==False:
            easyparam = (self.inputs.nx, self.inputs.ny,
                         self.inputs.tra_lon, self.inputs.tra_lat,
                         self.inputs.size_x, self.inputs.size_y,
                         self.inputs.rot)
            self.display(''.join(('--- saving grid (nx=%i, ny=%i, ',
                                  'lon=%.2f, lat=%.2f, sx=%.2f, sy=%.2f, ',
                                  'rot=%.2f)\n')) % easyparam)
        else:
            easyparam = (self.outputs.h.shape[1]-2, self.outputs.h.shape[0]-2,
                         self.prt_grd[2], 
                         self.prt_grd[3],self.prt_grd[4],self.prt_grd[5],self.prt_grd[6])
            self.display(''.join(('--- saving grid (nx=%i, ny=%i,',
                                  'coef=%i, imin=%i, imax=%i, jmin=%i, ',
                                  'jmax=%i)\n')) % easyparam)

        #self.display('--- saving grid (nx=%i, ny=%i, lon=%i, lat=%i, sx=%i, sy=%i, rot=%i)\n' % easyparam)
        self.display('Saving to netcdf done')

