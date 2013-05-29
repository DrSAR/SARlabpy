# -*- coding: utf-8 -*-
"""
Created on Wed May 15 14:30:17 2013

@author: fmoosvi
"""

import sarpy
import sarpy.fmoosvi.analysis
import numpy
import pylab
import os
import json

def get_bbox(value,data_label,type=None):
    
    data = sarpy.Scan(value[data_label][0])
    shape = data.pdata[0].data.shape
    
    bbox = numpy.array([float(x) for x in value[data_label][1]])    
    bbox_px = (bbox.reshape(2,2).T*shape[0:1]).T.flatten()
    
    #TODO: this will need to be updated for python 3.x+
    bbox_px = map(int,bbox_px) # Casts all elements to ints
    
    if type is None:
        return bbox_px
    
    elif type == 'pct':
        return bbox
        
        
        

def generate_VTC(masterlist_name, data_label):
  
    mdata = os.path.expanduser(os.path.join('~','mdata',masterlist_name + '.json'))
    
    with open(mdata,'r') as master_file:
        master_list = json.load(master_file)
    
    for k,v in master_list.iteritems():
        
        try:
            data = sarpy.Scan(v[data_label][0])
            bbox_px = get_bbox(v,data_label)
    
            # Normalize data
            ndata = sarpy.fmoosvi.analysis.h_normalize_dce(data)

            # Set bounding boxes and get ready to join together
            ndata[bbox_px[0]:bbox_px[1],bbox_px[2]:bbox_px[3],:,:] = numpy.nan
            ndata[:,:,:,-1] = numpy.nan
    
            # Get useful params        
            x_size = ndata.shape[0]
            y_size = ndata.shape[1]
            num_slices = ndata.shape[2]
            reps = ndata.shape[3]
    
            # Reshape it  to stitch together all the data
            nrdata = numpy.empty([x_size,y_size*reps,num_slices])
            
            for s in xrange(num_slices):
                nrdata[:,:,s] = ndata[:,:,s,:].reshape([x_size,y_size*reps])
                
            data.store_adata(key='VTC', data = nrdata)
            
        except IOError:
            
            pass   


masterlist_name = 'HerP2'
data_label = '24h-DCE'

a = generate_VTC(masterlist_name, data_label)

def plotVTC(value, data_label, adata_label = None):

    mdata = os.path.expanduser(os.path.join('~','mdata',masterlist_name + '.json'))
    
    data = sarpy.Scan(value[data_label][0])
    bbox_px = get_bbox(value,data_label)
    
    if adata_label is None:
        
        dffd
        
    else:
        
        nrdata = sarpy.Scan(value[adata_label][0])
 
    x_bbox_px = bbox_px[1] - bbox_px[0] +1

    # Handle the special case of a corner bbox at origin
    if bbox_px[0]==0: x_bbox_px - 1
    
#    pylab.figure(figsize = [10,data.method.PVM_FovCm[1]/data.method.PVM_FovCm[0]*10])
#    pylab.imshow(data.pdata[0].data[:,:,3,80])
#    
    G = pylab.matplotlib.gridspec.GridSpec(x_bbox_px, 1)
    fig = pylab.figure(figsize = [10,data.method.PVM_FovCm[1]/data.method.PVM_FovCm[0]*10])       

    i = 0
    
    for x in xrange(bbox_px[0],bbox_px[1]):
     
        fig.add_subplot(G[i],frameon=False,xticks=[],yticks=[])
    
        d = nrdata[x,:,3]
    
        pylab.plot(d,'-',linewidth=0.3)
        pylab.ylim([-0.2,2])
        
        i+=1 # increment counter for the subplot (x)     




























#bounding_box,box_size = create_bounding_box(data,top_left)
#
#G = pylab.matplotlib.gridspec.GridSpec(box_size[0], box_size[1])
#i = 0
#
#for x in xrange(bounding_box[1],bounding_box[3]):
#    
#    j=0 # Reset counter for next x line
#    
#    for y in range(bounding_box[0],bounding_box[2]):
#        
#        fig.add_subplot(G[j,i],frameon=False,xticks=[],yticks=[])
#        
#        d = data[y,x,3,:]
#        
#        pylab.plot(d,'-',linewidth=0.3)
#        pylab.ylim([-0.2,2])
#        j+=1 # increment counter for the subplot (y) 
#    
#    i+=1 # increment counter for the subplot (x) 
#        
#fig.savefig('youwinfiras.pdf')
##
#
#
#


#Nbig = 6
#Nsmall = 3
#big = np.arange(36).reshape([6,6])
#array([[ 0,  1,  2,  3,  4,  5],
#       [ 6,  7,  8,  9, 10, 11],
#       [12, 13, 14, 15, 16, 17],
#       [18, 19, 20, 21, 22, 23],
#       [24, 25, 26, 27, 28, 29],
#       [30, 31, 32, 33, 34, 35]])
#
#small = big.reshape([Nsmall, Nbig/Nsmall, Nsmall, Nbig/Nsmall]).mean(3).mean(1)
#
#array([[  3.5,   5.5,   7.5],
#       [ 15.5,  17.5,  19.5],
#       [ 27.5,  29.5,  31.5]])


#def rebin(a, new_shape):
#    """
#    Resizes a 2d array by averaging or repeating elements, 
#    new dimensions must be integral factors of original dimensions
#    
#     Credit: https://gist.github.com/zonca/1348792
# 
#    Parameters
#    ----------
#    a : array_like
#        Input array.
#    new_shape : tuple of int
#        Shape of the output array
# 
#    Returns
#    -------
#    rebinned_array : ndarray
#        If the new shape is smaller of the input array, the data are averaged, 
#        if the new shape is bigger array elements are repeated
# 
#    See Also
#    --------
#    resize : Return a new array with the specified shape.
# 
#    Examples
#    --------
#    >>> a = np.array([[0, 1], [2, 3]])
#    >>> b = rebin(a, (4, 6)) #upsize
#    >>> b
#    array([[0, 0, 0, 1, 1, 1],
#           [0, 0, 0, 1, 1, 1],
#           [2, 2, 2, 3, 3, 3],
#           [2, 2, 2, 3, 3, 3]])
# 
#    >>> c = rebin(b, (2, 3)) #downsize
#    >>> c
#    array([[ 0. ,  0.5,  1. ],
#           [ 2. ,  2.5,  3. ]])
# 
#    """
#    m, n = new_shape    
#    
#    
#    if size(a.shape) == 2: 
#        M, N = a.shape
#        if m<M:
#            return a.reshape((m,M/m,n,N/n)).mean(3).mean(1)
#        else:
#            return numpy.repeat(numpy.repeat(a, m/M, axis=0), n/N, axis=1)    
#
#    elif size(a.shape) == 3: 
#        M, N, O = a.shape
#        o = O
#        if m<M:
#            return a.reshape((m,M/m,n,N/n,o,O/o)).squeeze().mean(3).mean(1)
#        else:
#            return numpy.repeat(numpy.repeat(a, m/M, axis=0), n/N, o/O, axis=1)    
#
#    elif size(a.shape) == 4: 
#        M, N, O, P = a.shape
#        o = O
#        p = P
#        if m<M:
#            return a.reshape((m,M/m,n,N/n,o,O/o, p , P/p)).squeeze().mean(3).mean(1)
#        else:
#            return numpy.repeat(numpy.repeat(a, m/M, axis=0), n/N, o/O, p/P, axis=1)     
 
 
#if __name__ == '__main__':
#    import doctest
#    doctest.testmod()
    
    
#bata = rebin(data[:,:,:,:], (16,32))

#   num_slices = data_list[0].shape[-1]
#    data_num = len(data_list)
#    
#    fig = pylab.figure(figsize = (14,3))
#    G = pylab.matplotlib.gridspec.GridSpec(data_num,num_slices)   
#    #clims = sarpy.fmoosvi.getters.get_image_clims(data_list[0])
#    
#    for slice in xrange(num_slices):
#        
#        for n in xrange(data_num):    
#            fig.add_subplot(G[n,slice],frameon=False, xticks=[], yticks=[])
#            
#            try:
#                data = data_list[n][0,:,:,slice]  
#            except:
#                data = data_list[n][:,:,slice] 
#            
#            a = pylab.imshow(data, cmap = colour_map )
#            
#            if isinstance(clims[0],(int, float, long)):                
#                a.set_clim(clims)
#            else:
#                a.set_clim(600,3500)
#            if n == 0:
#                pylab.title('Slice {0}'.format(slice+1), fontsize = 14)
#    
#
#    # Colobar set up
#    cax = fig.add_axes([0.9, 0.10, 0.01, 0.7])
#    cax.set_title(key_list[1], fontsize = 12)       
#    fig.colorbar(a, cax=cax)
#    
#    # Figure spacing adjustments
#    #fig.subplots_adjust(right = 0.85, wspace = 0.0001, hspace=0.0001)
#    G.tight_layout(fig, h_pad = 0.1, w_pad = 0.001)
#    G.update(right = 0.87)