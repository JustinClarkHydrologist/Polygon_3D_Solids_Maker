# Polygon_3D_Solids_Maker

Generate 3D solids from 2D polygon GIS features.

You can use this lightweight Python tool to convert 2D polygons into 3D objects using specified top and bottom elevations.

Designed to work with shapefiles and geopackages.

Tool supports export formats compatible with both **ArcGIS Pro** (as multipatch features) and **QGIS** (via 3D model formats like DAE or GLTF).

## Purpose

This tool helps environmental scientists, hydrologists, geologists, and GIS users visualize map-based features in 3D 
by extruding 2D polygons into vertical surfaces that represent structural or hydrogeologic facets.

## Features

- Read input polygons from **shapefiles** or **geopackages**
- Extrude to create vertical walls using a user-defined top and bottom elevation
- Export to:
  - **Multipatch** (ArcGIS Pro)
  - **DAE / GLTF** (QGIS & 3D modeling software)

> Future Update: Support for variable vertex-wise extrusion based on elevation data from DEMs or point clouds.

## Requirements

- Python 3.9+
- [GeoPandas](https://geopandas.org/)
- [Shapely](https://shapely.readthedocs.io/)
- [Trimesh](https://trimsh.org/)
- [Fiona](https://fiona.readthedocs.io/)
- (Optional) `arcpy` for exporting to multipatch (ArcGIS Pro only)
- (Optional) `pycollada`, `pyvista`, or `vtk` for QGIS-compatible 3D formats

## Usage

```python

# Coming soon: Script-based CLI and example input/output files

Example file from Texas Water Development Board - GAM Downloads:
https://gw-models.s3.amazonaws.com/Download_GAMs/glfc_n/glfc_n_v4.01_GULF2023_modelGrid.7z
GULF_modelGrid.gdb

