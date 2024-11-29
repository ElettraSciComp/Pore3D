
import pypore3d.p3d_SITK_read_raw
from pypore3d.p3d_SITK_read_raw import *

import pypore3d.p3dFiltPy
from pypore3d.p3dFiltPy import py_p3dReadRaw8,py_p3dWriteRaw8


###### common 8 bit
def sitk_to_p3d_file_format (sitk_image, dimx, dimy, dimz):
    sitk.WriteImage(sitk_image , "tempMeta.mhd")
    im_P3D_format_out = py_p3dReadRaw8( 'tempMeta.raw', dimx,dimy,dimz)
    os.remove('tempMeta.raw')
    return im_P3D_format_out

def p3d_to_sitk_file_format (p3d_image, dimx,dimy,dimz):
    py_p3dWriteRaw8( p3d_image, 'tempOutMeta.raw',dimx,dimy,dimz)
    im_SITK_format_out = Read_Raw8("tempOutMeta.raw", [dimx,dimy,dimz])
    os.remove('tempOutMeta.raw')
    return im_SITK_format_out

def apply_rescaler(filtered_img,dimx,dimy,dimz):
    rescaler = sitk.RescaleIntensityImageFilter()
    rescaler.SetOutputMinimum(0)
    rescaler.SetOutputMaximum(255)
    rescaledImg = rescaler.Execute(filtered_img)
    return rescaledImg
