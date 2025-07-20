"""     Python based 3D object creator
             Started  on  4/22/2025
             Last Updated 7/20/2025 (July 20th, I am from the US)

@author: Justin A. Clark

This program takes a data from a GIS polygon feature and converts that feature into a 3D object
   - Can read shapefile and geopackages
   - Load MODFLOW Model Grids
   - GIS output formats will be selected by User
   
Test Run Time = X min XX Seconds

Python version: 3.9.13
"""
import geopandas as gpd                   # For reading shapefiles or geopackages
from shapely.geometry import Polygon      # To check if a feature is a polygon
import trimesh                            # For creating and exporting 3D mesh geometry
import os                                 # For file and directory operations
from pathlib import Path                  # For more readable and safer file path management

## Find kernel's current directory, update to desired location if necessary
## Windows and Linux path are formatted differently
path = r'/home/jac/Mdl/TX/GulfCoast_North/ModelGrid/GULF_modelGrid_QGIS.gdb'
#path = r'C:\MyPy\programs\Polygon_3D_Solids_Maker\Test01_Square_Phoenix.shp'  ##This is the test file I have provided for simple a 2D polygon shapefile
## TODO Implement Args parser or tkinter to get file path and name

# Function to preprocess grid data
def preprocess_geodata(input_file):
    """
    Load the input polygon data and determine the number of layers
    by parsing the column names. Supports both shapefiles and geopackages.

    Returns:
        gdf (GeoDataFrame): loaded GeoDataFrame
        top_field (str): name of the top elevation field for layer 1
        bot_field (str): name of the bottom elevation field for the last layer
        Layers (int): total number of model layers
    """
    # Detect file type and load accordingly
    if input_file.endswith(".shp") or input_file.endswith(".gpkg") or input_file.endswith(".gdb"):
        gdf = gpd.read_file(input_file)
    else:
        raise ValueError("Unsupported file format. Currenly only supports these file types: .shp, .gpkg, .gdb.")

    # List of all columns
    col_names = list(gdf.columns)

    #Example col_names: 
"""
[OBJECTID	row	col	natlRow	natlCol	top_1	bot_1	bot_2	bot_3	bot_4	bot_5	bot_6	IB_1	IB_2	IB_3	IB_4	IB_5	IB_6
  ZB_1	ZB_2	ZB_3	ZB_4	ZB_5	ZB_6	GHB	DRN_unname	DRN_named	DRN_Head	RIV_named	NHD_length	RIV_stage	RIV_bot
  fine_m2	fine_m3	fine_m4	fine_m5	rnb2	rnb3	rnb4	rnb5	IC_1	GCD	GMAnum	CNTY_NM	RWPA	X	Y	IB_active	Shape_Length	Shape_Area]
"""
    # Acceptable field name patterns for top of layer 1
    top_candidates = [c for c in col_names if c.lower() in ["top_1", "top1"]]

    # Raise error if top not found
    if not top_candidates:
        raise ValueError("Field for top of layer 1 not found. Acceptable options: top_1, Top_1, top1, Top1.")

    # Use the first match as the top field name
    top_field = top_candidates[0]

    # Find the bottom layer field using prefix logic (e.g., bot_1, bot2, BOT_10).
    bot_fields = []
    for col in col_names:
        col_lower = col.lower()
        if col_lower.startswith("bot") or col_lower.startswith("bot_"):
            # Extract numeric suffix
            suffix = ''.join(filter(str.isdigit, col))
            if suffix.isdigit():
                bot_fields.append((int(suffix), col))

    # === Identify bottom-most layer field ===
    # We previously created a list of tuples: (layer_number, field_name)
    # Now, we'll iterate through it to find the tuple with the highest layer number

    max_bot_layer_num = -1      # Initialize with invalid number
    bot_field = None            # Placeholder for field name with deepest layer

    for layer_num, field_name in bot_fields:
        if layer_num > max_bot_layer_num:
            max_bot_layer_num = layer_num
            bot_field = field_name  # Keep the name associated with the largest number

    # Final total number of layers (assuming model starts at layer 1)
    Layers = max_bot_layer_num

    # Generate error if field is not found
    if not bot_field:
        raise ValueError("Unable to determine deepest bottom layer field.")

    print(f"Top field: {top_field}, Bottom field: {bot_field}, Total layers: {Layers}")

    return gdf, top_field, bot_field, Layers
   
# Function to extrude a 2D polygon vertically into a 3D object
def extrude_polygon_to_3d(polygon: Polygon, z_min: float, z_max: float) -> trimesh.Trimesh:
    """
    Extrudes a shapely 2D polygon to a 3D wall between specified base and top elevations.
    
    Parameters:
        polygon (Polygon): The input 2D polygon geometry
        z_min (float): Bottom elevation of the extrusion
        z_max (float): Top elevation of the extrusion
    
    Returns:
        trimesh.Trimesh: A 3D mesh object representing the extruded polygon
    """
    # Extract exterior coordinates from the polygon
    exterior_coords = list(polygon.exterior.coords)

    # Create 2D path from coordinates using line segments
    path_2d = trimesh.path.Path2D(
        entities=[trimesh.path.entities.Line([i, i + 1]) for i in range(len(exterior_coords) - 1)],  # Lines between consecutive points
        vertices=exterior_coords  # Vertices of the polygon
    )

    # Extrude the 2D path vertically into a 3D shape with height = z_max - z_min
    mesh = path_2d.extrude(height=z_max - z_min)

    # Shift the entire mesh vertically so it starts at z_min
    mesh.apply_translation([0, 0, z_min])

    return mesh  # Return the 3D mesh

# Function to process all polygons in the input file
def process_geodata(gdf, top_field: str, bot_field: str, output_dir: str):
    """
    Extrude each polygon from its layer-1 top to deepest bottom and
    write an STL mesh.  All elevations come from attribute fields.
    
    Parameters:
        input_file (str): Path to the shapefile or geopackage (.shp or .gpkg)
        top_z (float): Top elevation of extrusion
        bottom_z (float): Bottom elevation of extrusion
        output_dir (str): Directory to save the exported 3D models
    """
    os.makedirs(output_dir, exist_ok=True)

    for idx, row in gdf.iterrows():
        geom = row.geometry
        if isinstance(geom, Polygon):
            z_top = row[top_field]        # e.g. value from top_1
            z_bot = row[bot_field]        # e.g. value from bot_6

            # Skip rows with missing elevation data
            if z_top is None or z_bot is None:
                print(f"  -- Missing z values at index {idx}; skipped.")
                continue

            # Extrude polygon into 3D mesh
            mesh = extrude_polygon_to_3d(geom, z_min=z_bot, z_max=z_top)

            # Define output path for this mesh using index
            out_path = Path(output_dir) / f"layer1_full_{idx}.stl"

            # Export the mesh as an STL file
            mesh.export(str(out_path))

            # Print confirmation message
            print(f"Exported: {out_path}")

        else:
            # Skip non-polygon features (e.g., points, lines), print message indicating occurrence
            print(f"  -- Non-polygon geometry at index {idx}; skipped.")


    # Loop through each row (feature) in the GeoDataFrame
    for idx, row in gdf.iterrows():
        geom = row.geometry  # Extract geometry from the current row

        # Check if geometry is a valid Polygon
        if isinstance(geom, Polygon):
            # Extrude polygon into 3D mesh
            mesh = extrude_polygon_to_3d(geom, z_min=bottom_z, z_max=top_z)

            # Define output path for this mesh using index
            out_path = Path(output_dir) / f"extruded_polygon_{idx}.stl"

            # Export the mesh as an STL file
            mesh.export(str(out_path))

            # Print confirmation message
            print(f"Exported: {out_path}")

        else:
            # Skip non-polygon features (e.g., points, lines)
            print(f"Skipping non-polygon geometry at index {idx}")

# Main execution block for running this script directly
if __name__ == "__main__":
    input_file = path                              # use your path variable
    gdf, top_field, bot_field, Layers = preprocess_geodata(input_file)

    # Define the output folder to save 3D models
    output_folder = "outputs_layer1_full"          # change if you like
    print(f"\n>>> Extruding polygons from {top_field} down to {bot_field}")

    # Run the full processing function with specified inputs
    process_geodata(gdf, top_field, bot_field, output_folder)

##############################################################################
#  ### ### ### ### ### #### ### ### ### ### ### #### ### ### ### ### ### ###  #
##   Example and Test Code Used for This Program   ##
"""




"""

###############################################################################
#  ### ### ### ### ### #### ### ### ### ### ### #### ### ### ### ### ### ###  #
##   Websites Visited During Code Making  ##
"""

https://gis.stackexchange.com/questions/349309/visualize-polygonz-shapefiles-in-qgis


"""
