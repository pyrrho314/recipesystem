try:
    import matplotlib.pyplot as plt
    from scipy.ndimage.interpolation import zoom
    def mpl_img(innd, zoom_factor = None, out_shape = None):
        if out_shape:
            w,h = innd.shape
            ow,oh = out_shape
            zf_1 = float(ow)/w
            zf_2 = float(oh)/h 
            zf = min(zf_1, zf_2)
            nd = zoom(innd, zf, prefilter = False)
        elif zoom_factor:
            nd = zoom(innd, zoom_factor, prefilter=False)
        else:
            nd = innd
        a = plt.imshow(nd)
        return a
except:
    def mpl_img(*args, **argv):
        return None


#### older

try:
    import numpy as np
    from scipy.ndimage.interpolation import zoom
    from osgeo import gdal, gdalnumeric
    def gdal2png(nd,output_name,proto_image, zoom_factor= None):
        #nd.resize((200,200))
        if zoom_factor:
            nd = zoom(nd, zoom_factor, prefilter=False)
        
        nd = nd.astype(float)
        input_array = (nd - nd.min())/np.ptp(nd)*256
        
        colorrdc = np.zeros((3,input_array.shape[0],input_array.shape[1]),dtype=np.float16)
        x=0
        y=0
        #ndmax = input_array.max()
        #ndmin = input_array.min()
        #ndtop = ndmax - ndmin
        #input_array = (input_array-ndmin)/ndtop
        if False:
            for x in range(input_array.shape[0]):
                for y in range(input_array.shape[1]):
                    if input_array[x,y] < 0.5:
                        colorrdc[0,x,y] = input_array[x,y]
                        colorrdc[1,x,y] = input_array[x,y]
                        colorrdc[2,x,y] = 0
                    else:
                        colorrdc[0,x,y] = input_array[x,y]
                        colorrdc[1,x,y] = input_array[x,y]
                        colorrdc[2,x,y] = 0
            rescaled_rdc = colorrdc * 255
            rescaled_rdc = rescaled_rdc.astype(gdalnumeric.uint8)
        colorrdc[0] = input_array
        colorrdc[1] = input_array
        #colorrdc[2] = input_array
        rescaled_rdc = colorrdc.astype(gdalnumeric.uint8)
        gdalnumeric.SaveArray(rescaled_rdc, output_name, format="PNG",prototype=proto_image)
        return rescaled_rdc
except:
    def gdal2png(a,b,c):
        return None

try:
    from IPython.display import Image
    
    def ipython_load_img (filename):
        im = Image("test.png")
        return im
except:
    def ipython_load_img(filename):
        return None