import matplotlib
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.widgets import PolygonSelector
import numpy
import pickle
import ipywidgets
import ipywidgets.widgets
from IPython.display import display

def distance_sqrd_P2segment(P, P0, P1):
    '''distance squared of P from finite line segment P0->P1

	a little tricky since there are three cases to consider
    implementation http://geomalgorithms.com/a02-_lines.html'''
    v = P1 - P0 # line segment
    w = P - P0  # point to first point
    c1 = numpy.dot(w, v)
    if c1 <= 0: # proj outside and 'before' P0
        return numpy.sum((P-P0)**2)
    c2 = numpy.dot(v,v)
    if c2 <= c1: # proj outside and 'after' P1
        return numpy.sum((P-P1)**2)
    # foot of projection between P0 and P1 -> get length of projection
    b = c1 / c2
    Pb = P0 + b*v
    return numpy.sum((P-Pb)**2)

def DrawROIMultiSclice(scn, roilabel_start='polygon_default'):
    if not roilabel_start.startswith('polygon_'):
        raise ValueError('roilabels need to start with "polygon_" by convention')
    data = scn.pdata[0].data
    fig, ax = plt.subplots()
    nr_slices = data.shape[2]
    imageref = ax.imshow(data[:,:,nr_slices//2])

    class MyPolygonSelector(matplotlib.widgets.PolygonSelector):
        def __init__(self, ax, fhandle, max_ds=10):
            super().__init__(ax, fhandle)
            self.__showverts = True
            self.__max_ds = max_ds

        def get_ind_under_cursor(self, event):
            'get the index of the vertex under cursor if within max_ds tolerance'
            # display coords
            xy = numpy.asarray(lasso.verts)
            xt, yt = xy[:, 0], xy[:, 1]
            eventxy = lasso.ax.transData.inverted().transform((event.x, event.y))
            d = numpy.sqrt((xt - eventxy[0])**2 + (yt - eventxy[1])**2)
            indseq = numpy.nonzero(numpy.equal(d, numpy.amin(d)))[0]
            ind = indseq[0]
            if d[ind] >= self._MyPolygonSelector__max_ds:
                ind = None
            return ind

        def get_segment_under_cursor(self, event):
            '''get the index of the segment closest to cursor coordinates
            
            iteration over all segments to find the segment closest to clicked position'''
            # display coords -> data coordinates
            eventxy = lasso.ax.transData.inverted().transform((event.x, event.y))
            # lasso coords for n vertices
            P0 = numpy.asarray(lasso.verts)  # starting point
            P1 = numpy.roll(P0, -1, axis=0)  # end point of segments
            # loop over all segments to find distance
            dist = numpy.zeros_like(P0[:,0])
            for i in range(len(lasso.verts)):
                dist[i]=distance_sqrd_P2segment(eventxy, P0[i,:], P1[i,:])
            indseg = numpy.nonzero(numpy.equal(dist, numpy.amin(dist)))[0][0]
            return indseg, eventxy

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
        fig.canvas.draw_idle()
        store_ROI(verts_2_save)

    def store_ROI(verts_2_save):
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

    def press(event):
        '''handler of keypresses in lasso widget'''
        if event.key == '+':
            SliceWidget.value += 1
        elif event.key == '-':
            SliceWidget.value -= 1
        elif event.key =='t':
            lasso._MyPolygonSelector__showverts = not(lasso._MyPolygonSelector__showverts)
            lasso.line.set_visible(lasso._MyPolygonSelector__showverts)
        elif event.key =='d':
            if len(lasso.verts) > 3:
                ind = lasso.get_ind_under_cursor(event)
                if ind is not None:
                    updated_verts = [tup for i,tup in enumerate(lasso.verts) if i!=ind]
                    updated_verts = numpy.append(updated_verts, [updated_verts[0]], axis=0)
                    lasso.setverts(updated_verts)
                    store_ROI(updated_verts)
        elif event.key =='i':
            inds, coords = lasso.get_segment_under_cursor(event)
            if inds is not None:
                updated_verts = numpy.append(lasso.verts, [lasso.verts[0]], axis=0)
                updated_verts = numpy.insert(updated_verts, inds+1, coords, axis=0)
                lasso.setverts(updated_verts)
                store_ROI(updated_verts)

    fig.canvas.mpl_connect('key_press_event', press)

    # get me all polygon adata:
    roi_labels = [i for i in scn.adata if i.startswith('polygon_')]
    if not roi_labels:
        roi_labels = [roilabel_start]
    lasso = MyPolygonSelector(ax, onselect)

    SliceWidget = ipywidgets.widgets.IntSlider(min=0, max=nr_slices-1,
                                               value = nr_slices//2,
                                               description = 'Slice #')
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
