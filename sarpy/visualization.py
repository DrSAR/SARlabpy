import matplotlib
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.widgets import PolygonSelector
import numpy
import pickle
import ipywidgets
import ipywidgets.widgets
from IPython.display import display

class CompoundDropdownWidgetClass(ipywidgets.widgets.Dropdown):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.TextInput = ipywidgets.widgets.Text(description='New ROI',
                                        continuous_update=False,
                                        value='',
                                        placeholder='Type ROI label to be created',)
        self.ROI2Delete=False
        self.DeleteButton = ipywidgets.widgets.Button(description='Delete',
                                            tooltip='Delete current ROI label')
        def on_DeleteButton_clicked(b):
            self.ROI2Delete = self.value # remember what to delete
            self.options = tuple([a for a in self.options if a !=self.value])
            if self.options:
                self.value = self.options[0]
            else:
                self.options = ['default']
            change = {'name': 'value',
                      'new': None,
                      'old': None,
                      'type': 'change'}
        self.DeleteButton.on_click(on_DeleteButton_clicked)
        self.AddButton = ipywidgets.widgets.Button(description='Create',
                                            tooltip='Create ROI with New Name')
        def on_AddButton_clicked(b):
            if self.TextInput.value and (self.TextInput.value not in self.options):
                self.options = self.options + (self.TextInput.value,)
                self.value = self.TextInput.value
        self.AddButton.on_click(on_AddButton_clicked)
        self.HorizBox = ipywidgets.widgets.HBox([self, 
                                                 self.DeleteButton, 
                                                 self.TextInput, 
                                                 self.AddButton])

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

def DrawROIMultiSclice(scn, roilabel_start=None):
    data = scn.pdata[0].data
    fig, ax = plt.subplots()
    nr_slices = data.shape[2]
    imageref = ax.imshow(data[:,:,nr_slices//2])

    class MyPolygonSelector(matplotlib.widgets.PolygonSelector):
        def __init__(self, ax, fhandle, max_ds=10):
            super().__init__(ax, fhandle)
            self.__showverts = True
            self.__max_ds = max_ds
            self.__showing_mask = False

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

    x, y = numpy.meshgrid(numpy.arange(data.shape[1]), numpy.arange(data.shape[0]))
    pix = numpy.vstack((x.flatten(), y.flatten())).T

    def onselect(verts):
        # Select elements in original array bounded by selector path:
        p = Path(verts)
        # append first point to end to close polygon
        verts_2_save = numpy.append(verts, [verts[0]], axis=0)
        ind = p.contains_points(pix, radius=0)
        selected = numpy.zeros_like(data[:,:,0])
        selected.flat[ind] = numpy.max(data)
        lasso.__mask = selected
        fig.canvas.draw_idle()
        store_ROI(verts_2_save)

    def store_ROI(verts_2_save):
        roi_label = imageref.get_figure().canvas.__CurrentRoiLabel
        currentslice = imageref.get_figure().canvas.__CurrentSlice
        adata = scn.adata.get('polygon_'+roi_label)
        if adata is None:
		# create new dictionary
            vert_data = {currentslice: verts_2_save}
        else:
            vert_data = adata.data
            vert_data[currentslice] = verts_2_save
        scn.store_adata(key='polygon_'+roi_label, force=True, data=vert_data)
        
    def ROIplotUpdater(change):
        # there appear to be a lot of events on the widgets not all of which
        # warrant a response here ...
        curslice = SliceWidget.value
        # remember which slice we're on now
        imageref.get_figure().canvas.__CurrentSlice = curslice
        imageref.set_data(data[:,:,curslice])
        imageref.set_clim(vmin=WindowingWidget.value[0], 
                          vmax=WindowingWidget.value[1])
        imageref.set_cmap(CmapWidget.value)
        # remember which ROI we're on now
        if not RoiLabelWidget.options:
            RoiLabelWidget.options=['default']
        roi_label = RoiLabelWidget.value
        imageref.get_figure().canvas.__CurrentRoiLabel = roi_label
        # load vertices for current slice
        vert_data = scn.adata.get('polygon_'+roi_label) # return None if not there
        # show vertices
        if (vert_data is None) or (curslice not in vert_data.data):
            lasso.setverts(numpy.array([[], []]).T)
        else:
            lasso.setverts(vert_data.data[curslice])
        # is there a ROI that needs deleting?
        if RoiLabelWidget.ROI2Delete:
            scn.rm_adata('polygon_'+RoiLabelWidget.ROI2Delete)
            RoiLabelWidget.ROI2Delete= False

    def press(event):
        '''handler of keypresses in lasso widget'''
        if event.key == '+':
            SliceWidget.value += 1
        elif event.key == '-':
            SliceWidget.value -= 1
        elif event.key =='t':
            lasso._MyPolygonSelector__showverts = not(lasso._MyPolygonSelector__showverts)
            lasso.artists[1].set_visible(lasso._MyPolygonSelector__showverts)
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
        elif event.key == 'm':
            # make vertices and markers (in)visible
            lasso._MyPolygonSelector__showing_mask = not(lasso._MyPolygonSelector__showing_mask)
            for a in lasso.artists:
                    a.set_visible(not(lasso._MyPolygonSelector__showing_mask))
            lasso._MyPolygonSelector__showverts = not(lasso._MyPolygonSelector__showing_mask)
            if lasso._MyPolygonSelector__showing_mask:
                imageref.set_data(lasso.__mask)
            else:
                imageref.set_data(data[:,:,imageref.get_figure().canvas.__CurrentSlice])

    fig.canvas.mpl_connect('key_press_event', press)
    lasso = MyPolygonSelector(ax, onselect)

    SliceWidget = ipywidgets.widgets.IntSlider(min=0, max=nr_slices-1,
                                               value = nr_slices//2,
                                               description = 'Slice #')
    SliceWidget.observe(ROIplotUpdater, 'value')

    # get me all polygon adata and remove 'polygon_' from label:
    roilabels = [i[len('polygon_'):] for i in scn.adata if i.startswith('polygon_')]
    if roilabel_start is None and not(roilabels):
        roilabels = ['default']
    if roilabel_start in roilabels:
        roi_firstchoice = roilabel_start
    elif roilabel_start is None:
        roi_firstchoice = roilabels[0]
    else: #need to add the supplied label to the choices
        roilabels.append(roilabel_start)
    RoiLabelWidget = CompoundDropdownWidgetClass(options = roilabels,
                                                 value = roi_firstchoice,
                                                 description = 'Current ROI')
    RoiLabelWidget.observe(ROIplotUpdater,'value')

    WindowingWidget = ipywidgets.widgets.FloatRangeSlider(
        value=[numpy.percentile(data, 5), numpy.percentile(data, 95)],
        min=data.min(),
        max=data.max(),
        description='windowing:',
        continuous_update=True,
        orientation='horizontal',
        readout=True,
        readout_format='.1f',
        layout=ipywidgets.Layout(width='100%')
    )
    WindowingWidget.observe(ROIplotUpdater,'value')

    CmapWidget = ipywidgets.widgets.Dropdown(
                    options = ['Greys', 'gray', 'viridis', 'plasma', 'PiYG'],
                    value = 'gray',
                    description = 'Colormap')
    CmapWidget.observe(ROIplotUpdater,'value')

    my_combo_widget = ipywidgets.widgets.VBox([SliceWidget,
                                              RoiLabelWidget.HorizBox,
                                              WindowingWidget,
                                              CmapWidget])
    # initial call to handler so that widget defaults and display are in sync
    ROIplotUpdater(None)
    display(my_combo_widget)

