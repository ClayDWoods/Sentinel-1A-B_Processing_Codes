import rasterio
import os
import numpy as np
from rasterio import Affine
from rasterio.warp import reproject, Resampling
from rasterio.windows import get_data_window


def tif_2_array(tif_file):
    tif_open = rasterio.open(tif_file)
    tif_array = tif_open.read(1)
    tif_open.close()
    return tif_array


def tif_2_array_multiband(tif_file, band):
    tif_open = rasterio.open(tif_file)
    tif_array = tif_open.read(band)
    tif_open.close()
    return tif_array


def tif_2_array_window(tif_file, window):
    tif_open = rasterio.open(tif_file)
    tif_array = tif_open.read(1, window=window)
    tif_open.close()
    return tif_array


def get_kwargs(tif_file):
    tif_open = rasterio.open(tif_file)
    kwargs_tif = tif_open.meta.copy()
    tif_open.close()
    return kwargs_tif


def write_reprojected_array_2_tif(projected_array, new_tif_name_full_path, kwargs_tif):
    if new_tif_name_full_path.exists():
        os.remove(new_tif_name_full_path)
    reprojected_tif = rasterio.open(new_tif_name_full_path, 'w', **kwargs_tif)
    reprojected_tif.write(projected_array.astype('float32'), indexes=1)
    reprojected_tif.close()


def get_tif_width(tif_file):
    tif_open = rasterio.open(tif_file)
    tif_width = tif_open.width
    tif_open.close()
    return tif_width


def get_tif_height(tif_file):
    tif_open = rasterio.open(tif_file)
    tif_height = tif_open.height
    tif_open.close()
    return tif_height


def get_tif_transform(tif_file):
    tif_open = rasterio.open(tif_file)
    tif_transform = tif_open.transform
    tif_open.close()
    return tif_transform


def get_tif_crs(tif_file):
    tif_open = rasterio.open(tif_file)
    tif_crs = tif_open.crs
    tif_open.close()
    return tif_crs


def make_affine(a, b, c, d, e, f):
    new_affine = Affine(a, b, c, d, e, f)
    return new_affine


def reproject_tif_array(tif_file, master_height, master_width, old_array, old_transform, base_crs, master_transform):
    tif_file_open = rasterio.open(tif_file)
    reprojected_array = np.zeros((master_height, master_width))
    reproject(
        old_array, reprojected_array,
        src_transform=old_transform,
        dst_transform=master_transform,
        src_crs=base_crs,
        dst_crs=base_crs,
        resampling=Resampling.nearest)
    tif_file_open.close()
    return reprojected_array


def reproject_tif_array_ztd(tif_file_base, tif_file_2_reproject):
    width_ = get_tif_width(tif_file_base)
    height_ = get_tif_height(tif_file_base)
    base_crs = get_tif_crs(tif_file_base)
    base_transform = get_tif_transform(tif_file_base)
    array_2_reproject = tif_2_array(tif_file_2_reproject)
    transform_2_reproject = get_tif_transform(tif_file_2_reproject)
    reprojected_array = np.zeros((height_, width_))
    reproject(
        array_2_reproject,
        reprojected_array,
        src_transform=transform_2_reproject,
        dst_transform=base_transform,
        src_crs=base_crs,
        dst_crs=base_crs,
        resampling=Resampling.cubic)
    return reprojected_array.astype('float32')


def lon_lat_2_col_row(tif_file, lon, lat):
    affine_transform = get_tif_transform(tif_file)
    a, c = affine_transform.a, affine_transform.c
    e, f = affine_transform.e, affine_transform.f
    col = int(abs((lon - c) / a))
    row = int(abs((lat - f) / e))
    return col, row


def col_row_2_lon_lat(tif_file, col, row):
    affine_transform = get_tif_transform(tif_file)
    lon, lat = affine_transform * (col, row)
    return lon, lat


def lon_lat_2_col_row_using_transform(transform, lons, lats):
    a, c = transform[0], transform[2]
    e, f = transform[4], transform[5]
    col = int((lons - c) / a)
    row = int((lats - f) / e)
    return col, row


def col_row_2_lon_lat_waffine(master_affine, col, row):
    lon, lat = master_affine * (col, row)
    return lon, lat


def get_data_window_f(tif_file):
    tif_open = rasterio.open(tif_file)
    window = get_data_window(tif_open.read(1, masked=True))
    tif_open.close()
    return window


def windowed_kwargs(tif_file, window):
    kwargs_old = get_kwargs(tif_file)
    src_transform = get_tif_transform(tif_file)
    kwargs_old.update({
        'height': window.height,
        'width': window.width,
        'transform': rasterio.windows.transform(window, src_transform)})
    return kwargs_old


def tif_bounds_ztd(tif_file):
    tif_open = rasterio.open(tif_file)
    tif_bounds = tif_open.bounds
    left, bottom, right, top = (round(tif_bounds[0], 1) - 0.1), (round(tif_bounds[1], 1) - 0.1), (
                round(tif_bounds[2], 1) + 0.1), (round(tif_bounds[3], 1) + 0.1)
    tif_open.close()
    return left, bottom, right, top


def return_bounds(tif_file):
    tif_open = rasterio.open(tif_file)
    tif_bounds = tif_open.bounds
    left, bottom, right, top = tif_bounds[0], tif_bounds[1], tif_bounds[2], tif_bounds[3]
    return left, bottom, right, top


def find_bounds_from_affine_w_h(affine, width, height):
    upper_left_lon, upper_left_lat = col_row_2_lon_lat_waffine(affine, 0, 0)
    lower_right_lon, lower_right_lat = col_row_2_lon_lat_waffine(affine, width-1, height-1)
    left_, right_, bot_, top_ = upper_left_lon, lower_right_lon, lower_right_lat, upper_left_lat
    return left_, bot_, right_, top_


def make_affine_from_list(transform_list):
    new_affine = Affine(transform_list[0], transform_list[1], transform_list[2], transform_list[3], transform_list[4], transform_list[5])
    return new_affine