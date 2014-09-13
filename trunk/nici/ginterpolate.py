from scipy import ndimage as nd
import numpy as np

def ginterpolate (im, X, Y=None, output_type=None, order=1, 
                  mode='constant', cval=0.0, q=False) :
    """
      Emulate the IDL interface to the intepolate
      function using the ndimage.map_coordinates call
      
      @type im:   2-D image buffer
      @type X:    2-D index array with the x-location for which 
                  interpolation is desired. 
      @type Y:    2-D index array with the y-location for which 
                  interpolation is desired. 
      @type order: Default:1. Order of interpolation. The default for
                   map_coordinates is 3
      @type q:     Default False. If True, there is only X and
                   we use only the upper left quarter of the symmetric X
                   (w/r to the array center). NOTE: Currently it is slower
                   to use this option. 

      For more information see help ndimage.map_coordinates

    """
    dim = np.size(np.shape(im))
    szx = np.size(np.shape(X))
    X = np.array(X)
    if Y != None: Y = np.array(Y)

    if dim == 1:
       if szx == 1:
           X = X.reshape([1,np.size(X)])
           z = nd.map_coordinates(im,X, output=output_type, order=order,
                     mode=mode, cval=cval)
       elif szx == 2:
           lim = np.size(im)
           limy,limx = lim,lim
           if Y == None: 
               Y = X
               coords = np.array([X,Y])
               if q:
                   nx,ny = np.shape(X)
                   nx,ny=nx/2,ny/2
                   coords = np.array([X[:ny,:nx],X[:ny,:nx]])
                   limy,limx = lim/2,lim
           else:
               coords = np.array([X,Y])
           z = nd.map_coordinates(np.resize(im,[limy,limx]),coords,
                     output=output_type, order=order, mode=mode, cval=cval)
    elif dim == 2:
       if szx == 1 and Y != None:
           z = nd.map_coordinates(im,[X,Y],output=output_type, order=order,
                     mode=mode, cval=cval)
       elif szx == 2:
           if Y == None: Y = X
           z = nd.map_coordinates(im,[X,Y],output=output_type, order=order,
                     mode=mode, cval=cval)
    else:
       print 'ERROR, dimension of input buffer not supported: ',len
       return 

    return z

    
