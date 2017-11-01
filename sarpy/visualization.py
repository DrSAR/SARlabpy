import matplotlib
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.widgets import PolygonSelector
import numpy
import pickle
import ipywidgets
import ipywidgets.widgets
from IPython.display import display

def DrawROIMultiSclice(scn, roilabel_start='polygon_default'):
    if not roilabel_start.startswith('polygon_'):
        raise ValueError('roilabels need to start with "polygon_" by convention')
    data = scn.pdata[0].data
    fig, ax = plt.subplots()
    nr_slices = data.shape[2]
    imageref = ax.imshow(data[:,:,nr_slices//2])

    class MyPolygonSelector(matplotlib.widgets.PolygonSelector):
        def setverts(self,verts):
            self._xs = [row[0] for row in verts]
            self._ys = [row[1] for row in verts]
            self._polygon_completed = True
            self._draw_polygon()

    def onselect(verts):
        # Select elements in original array bounded by selector path:
        p = Path(verts)
        # append first point to end to close polygon
        verts_2_save = numpy.append(verts, [verts[0]], axis=0)
        #ind = p.contains_points(pix, radius=10)
        fig.canvas.draw_idle()
        
        roi_label = imageref.get_figure().canvas.__CurrentRoiLabel
        currentslice = imageref.get_figure().canvas.__CurrentSlice
        adata = scn.adata.get(roi_label)
        if adata is None:
            # create new dictionary
            vert_data = {currentslice: verts_2_save}
        else:
            vert_data = adata.data
            vert_data[currentslice] = verts_2_save
        scn.store_adata(key=roi_label, force=True, data=vert_data)

    def f(n, roi_label, windowing, cmap):
        # update image display
        imageref.set_data(data[:,:,n])
        imageref.set_clim(vmin=windowing[0], vmax=windowing[1])
        imageref.set_cmap(cmap)
        # remember which slice we're on now
        imageref.get_figure().canvas.__CurrentSlice = n
        imageref.get_figure().canvas.__CurrentRoiLabel = roi_label
        # load vertices for current slice
        vert_data = scn.adata.get(roi_label) # return None if not there
        # show vertices
        if (vert_data is None) or (n not in vert_data.data):
            lasso.setverts(numpy.array([[], []]).T)
        else:
            lasso.setverts(vert_data.data[n])

    lasso = MyPolygonSelector(ax, onselect)
    SliceWidget = ipywidgets.widgets.IntSlider(min=0, max=nr_slices-1,
                                               value = nr_slices//2,
                                               description = 'Slice #')
    # get me all polygon adata:
    roi_labels = [i for i in scn.adata if i.startswith('polygon_')]
    if roilabel_start not in roi_labels:
        roi_labels.append(roilabel_start)
    RoiLabelWidget = ipywidgets.widgets.Dropdown(options = roi_labels, 
							value=roilabel_start, 
							description = 'ROI labels:')
    WindowingWidget = ipywidgets.widgets.FloatRangeSlider(
    		value=[numpy.percentile(data, 5), numpy.percentile(data, 95)],
    		min=data.min(), max=data.max(),
    		description='windowing:',
    		continuous_update=True,
    		orientation='horizontal',
    		readout=True,
    		layout=ipywidgets.Layout(width='90%'))
    CmapWidget = ipywidgets.widgets.Dropdown(
                options = ['Greys', 'gray', 'viridis', 'plasma', 'PiYG'],
                value = 'gray',
                description = 'Colormap')

    my_plts = ipywidgets.widgets.interactive(f,
                                         n=SliceWidget,
                                         roi_label=RoiLabelWidget,
                                         windowing=WindowingWidget,
                                         cmap=CmapWidget)
    display(my_plts)
