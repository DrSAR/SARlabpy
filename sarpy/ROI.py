'''Tools to manipulate region of interest data'''
import numpy
import matplotlib.path
try:
    import cv2
except ImportError:
    print('error importing OpenCV package cv2 -- this means PolygonFromMask will fail')

def MaskFromPolygon(verts, shape):
    x, y = numpy.meshgrid(numpy.arange(shape[1]), numpy.arange(shape[0]))
    pix = numpy.vstack((x.flatten(), y.flatten())).T
    p = matplotlib.path.Path(verts)
    ind = p.contains_points(pix, radius=-1)
    selected2D = numpy.zeros(shape=shape[0:2])
    selected2D.flat[ind] = 1
    return selected2D

def MaskFromADataPolygon(AData):
    '''convert multi-slice 2D polygons into binary masks'''
    shape = AData.parent.data.shape
    selected = numpy.zeros(shape=shape)
    for sliceidx, verts in AData.data.items():
        selected[:,:,sliceidx] = MaskFromPolygon(verts, shape[0:2])
    return selected

def PolygonFromMask(mask):
    im2, contours, hierarchy = cv2.findContours(mask.astype(numpy.uint8),
                                    cv2.RETR_TREE, 
                                    cv2.CHAIN_APPROX_TC89_L1,
#                                    cv2.CHAIN_APPROX_SIMPLE,
                                    )
    return contours[0].squeeze() 

