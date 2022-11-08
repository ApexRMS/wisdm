import rasterio
import rioxarray
# import numpy as np
from dask.distributed import Lock

def parc(inputFile, templateRaster, outputFile, chunks, resampling, mask, client):
    
    # set output structure 
    vrt_options = {
    'resampling': resampling,
    'crs': templateRaster.rio.crs ,
    'transform': templateRaster.rio.transfrom(),
    'height': templateRaster.rio.height,
    'width': templateRaster.rio.width,
    }
        
    # open and warp input raster
    with rasterio.open(inputFile) as src:
        with rasterio.vrt.WarpedVRT(src, **vrt_options) as vrt:
            ds = rioxarray.open_rasterio(vrt, chunks=chunks, lock=False)
            
            # clip to mask
            ds = ds.rio.clip(mask)

            # set NA values
            # if not np.issubdtype(ds.dtype, np.unsignedinteger):
            #     ds.rio.write_nodata(-9999, inplace = True)

            # write to disk
            ds.rio.to_raster(outputFile, tiled=True, lock=Lock("rio", client=client), windowed=True, overwrite=True, compress='lzw')

