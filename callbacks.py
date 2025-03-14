import base64
import io
import json
import re
import ssl
import certifi
import os
import requests
import urllib3
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
import plotly.io as pio
from io import StringIO
import dash
from dash.dependencies import Input, Output
from dash import Input, Output, State, dcc, html, dash_table
from dash.exceptions import PreventUpdate
import dash_bio as dashbio
from dash import ctx
import phylo_map
from Bio import Phylo, SeqIO
import plotly.colors as pcolors
from plotly.colors import qualitative
from dotenv import load_dotenv  

# Geopy (for geocoding city names)
import geopy.geocoders
from geopy.geocoders import Nominatim

# ✅ Force Python to use certifi's certificates
os.environ['SSL_CERT_FILE'] = certifi.where()
os.environ['SSL_CERT_DIR'] = certifi.where()

# ✅ Create an SSL context that explicitly uses certifi
ssl_context = ssl.create_default_context(cafile=certifi.where())

# ✅ Configure urllib3 to use the correct SSL certificates
http = urllib3.PoolManager(cert_reqs="CERT_REQUIRED", ca_certs=certifi.where())

# ✅ Suppress SSL warnings (prevents flooding logs with SSL errors)
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

load_dotenv()
print(f"Loaded API Key: {os.getenv('OPENCAGE_API_KEY')}")

# Store Markers in Memory
MARKERS = []
STANDALONE_MARKERS = []

def get_city_coordinates(city_name):
    """Fetch city coordinates using OpenCage API."""
    api_key = os.getenv("OPENCAGE_API_KEY")  # Load from environment variable
    if not api_key:
        return None, None, "⚠️ Missing OpenCage API Key."

    base_url = "https://api.opencagedata.com/geocode/v1/json"
    params = {"q": city_name, "key": api_key, "limit": 1}

    try:
        response = requests.get(base_url, params=params, verify=False, timeout=10)  # ✅ Disable SSL verification
        response.raise_for_status()

        data = response.json()
        if not data["results"]:
            return None, None, "⚠️ City not found. Please enter a valid city name."

        lat = data["results"][0]["geometry"]["lat"]
        lon = data["results"][0]["geometry"]["lng"]
        return float(lat), float(lon), None  # ✅ Success

    except requests.exceptions.RequestException as e:
        return None, None, f"⚠️ Error fetching city coordinates: {str(e)}"



def get_rectangular_coordinates(tree):
    xcoords = tree.depths(unit_branch_lengths=True)
    ycoords = {}

    def calc_y_coordinates(clade, current_y):
        if clade.is_terminal():
            ycoords[clade] = current_y
            return current_y + 1

        for subclade in clade:
            current_y = calc_y_coordinates(subclade, current_y)

        ycoords[clade] = (ycoords[clade.clades[0]] + ycoords[clade.clades[-1]]) / 2
        return current_y

    calc_y_coordinates(tree.root, 0)
    return xcoords, ycoords

def draw_clade_rectangular(clade, x_start, line_shapes, x_coords, y_coords):
    x_end = x_coords[clade]
    y_current = y_coords[clade]

    # Draw horizontal line for the branch
    line_shapes.append(dict(
        type='line',
        x0=x_start, y0=y_current, x1=x_end, y1=y_current,
        line=dict(color='black', width=1)
    ))

    # Draw vertical connecting lines for children
    if clade.clades:
        y_top = y_coords[clade.clades[0]]
        y_bottom = y_coords[clade.clades[-1]]
        line_shapes.append(dict(
            type='line',
            x0=x_end, y0=y_bottom, x1=x_end, y1=y_top,
            line=dict(color='black', width=1)
        ))

        for subclade in clade:
            draw_clade_rectangular(subclade, x_end, line_shapes, x_coords, y_coords)


def generate_location_colors(locations):
    unique_locations = locations.unique()
    colors = px.colors.qualitative.Plotly  # Pick a color palette
    color_map = {loc: colors[i % len(colors)] for i, loc in enumerate(unique_locations)}
    return color_map


def create_tree_plot(tree_file, metadata_file, show_tip_labels, height=None, width=None):
    """Generates a rectangular phylogenetic tree plot with MLST heatmap, bootstrap support, and legends with group titles."""

    # Load the tree
    tree = Phylo.read(tree_file, 'newick')
    tree.root_at_midpoint()

    # Load metadata
    metadata = pd.read_csv(metadata_file, sep='\t')
    if 'taxa' not in metadata.columns or 'location' not in metadata.columns or 'MLST' not in metadata.columns:
        raise ValueError("Metadata file must contain 'taxa', 'location', and 'MLST' columns.")

    metadata['location'] = metadata['location'].fillna('Unknown')
    metadata['MLST'] = metadata['MLST'].fillna('Unknown')

    # Generate location colors
    location_colors = {loc: px.colors.qualitative.Plotly[i % len(px.colors.qualitative.Plotly)]
                       for i, loc in enumerate(metadata['location'].unique())}

    # Generate MLST colors
    mlst_colors = {mlst: pcolors.qualitative.Vivid[i % len(pcolors.qualitative.Vivid)]
                   for i, mlst in enumerate(metadata['MLST'].unique())}

    # Compute tree node coordinates
    x_coords = tree.depths(unit_branch_lengths=True)
    y_coords = {}
    max_y = 0

    def assign_coordinates(clade, x_start=0, y_start=0):
        nonlocal max_y
        branch_length = clade.branch_length if clade.branch_length else 0.0
        x_current = x_start + branch_length

        if clade.is_terminal():
            x_coords[clade] = x_current
            y_coords[clade] = y_start
            max_y = max(max_y, y_start)
            return y_start + 1
        else:
            y_positions = []
            for child in clade.clades:
                y_start = assign_coordinates(child, x_current, y_start)
                y_positions.append(y_start - 1)

            x_coords[clade] = x_current
            y_coords[clade] = sum(y_positions) / len(y_positions)
            return y_start

    assign_coordinates(tree.root)

    # Set dynamic figure size based on the tree structure
    num_tips = sum(1 for _ in tree.get_terminals())  # Number of leaf nodes (taxa)
    max_label_length = max((len(clade.name) for clade in tree.get_terminals() if clade.name), default=10)

    # Dynamically adjust width and height
    height = max(800, num_tips * 25)  # Ensure enough space for tips
    width = max(1000, 800 + (max_label_length * 10))  # Ensure enough space for long labels

    # Set position for MLST heatmap squares
    mlst_x_position = max(x_coords.values()) + 0.02  

    # Lists for plot traces
    tree_line_traces = []
    bootstrap_markers = []
    tip_markers = []
    mlst_markers = []

    seen_locations = set()  
    seen_mlst = set()  

    # ** Add Legend Titles as Dummy Scatters **
    location_legend_title = go.Scatter(
        x=[None], y=[None], mode="markers",
        marker=dict(size=0, opacity=0),
        name="<b>Location</b>",
        showlegend=True
    )

    mlst_legend_title = go.Scatter(
        x=[None], y=[None], mode="markers",
        marker=dict(size=0, opacity=0),
        name="<b>MLST</b>",
        showlegend=True
    )

    # ** Draw Tree Branches as Scatter Lines **
    for clade in tree.find_clades(order='level'):
        x_start = x_coords[clade]
        y_start = y_coords[clade]

        if clade.clades:
            y_positions = [y_coords[child] for child in clade.clades]
            # Vertical line
            tree_line_traces.append(go.Scatter(
                x=[x_start, x_start], y=[min(y_positions), max(y_positions)],
                mode='lines', line=dict(color='black', width=2), showlegend=False
            ))
            # Horizontal lines
            for child in clade.clades:
                x_end = x_coords[child]
                y_end = y_coords[child]
                tree_line_traces.append(go.Scatter(
                    x=[x_start, x_end], y=[y_end, y_end],
                    mode='lines', line=dict(color='black', width=2), showlegend=False
                ))

        # Add bootstrap support markers (black diamonds) for nodes with confidence > 0.9
        if clade.confidence and clade.confidence > 0.9:
            bootstrap_markers.append(go.Scatter(
                x=[x_start], y=[y_start], mode='markers',
                marker=dict(size=12, color='black', symbol='diamond'),
                hoverinfo='text', text=f"Bootstrap: {clade.confidence}",
                showlegend=False
            ))

    # ** Draw Tip Markers Last (To Ensure They Appear on Top) **
    for clade in tree.get_terminals():
        x, y = x_coords[clade], y_coords[clade]
        meta_row = metadata[metadata['taxa'] == clade.name]

        if not meta_row.empty:
            location = meta_row['location'].iloc[0]
            mlst_value = meta_row['MLST'].iloc[0]
            color = location_colors.get(location, 'gray')
            mlst_color = mlst_colors.get(mlst_value, 'gray')

            show_location_legend = location not in seen_locations
            if show_location_legend:
                seen_locations.add(location)

            # **Conditionally render tip labels based on show_tip_labels**
            if show_tip_labels:
                tip_markers.append(go.Scatter(
                    x=[x], y=[y], mode='markers+text',
                    marker=dict(size=16, color=color, line=dict(width=2, color='black')),  
                    name=location if show_location_legend else None, 
                    text=f"{clade.name}",
                    textposition="middle right", textfont=dict(size=10), 
                    hoverinfo='text', showlegend=show_location_legend
                ))
            else:
                tip_markers.append(go.Scatter(
                    x=[x], y=[y], mode='markers',
                    marker=dict(size=16, color=color, line=dict(width=2, color='black')),  
                    name=location if show_location_legend else None, 
                    hoverinfo='text', showlegend=show_location_legend
                ))

            show_mlst_legend = mlst_value not in seen_mlst
            if show_mlst_legend:
                seen_mlst.add(mlst_value)

            mlst_markers.append(go.Scatter(
                x=[mlst_x_position], y=[y], mode='markers',
                marker=dict(size=20, color=mlst_color, symbol='square',
                    line=dict(width=2, color='black')),
                name=f"{mlst_value}" if show_mlst_legend else None,
                hoverinfo='text', text=f"MLST: {mlst_value}",
                showlegend=show_mlst_legend
            ))

    layout = go.Layout(
        title='Phylogenetic Tree with MLST Heatmap, Bootstrap Support, and Location Legend',
        xaxis=dict(title='Evolutionary Distance', showgrid=False, zeroline=False, range=[0, mlst_x_position + 0.01]),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-1, max_y + 1]),
        height=height, width=width, plot_bgcolor="rgb(240, 240, 250)"
    )

    # ** Order Matters: Draw Tree Lines First, Then Bootstrap, Then Tip Markers Last **
    return go.Figure(
        data=[location_legend_title] + tree_line_traces + bootstrap_markers + tip_markers + [mlst_legend_title] + mlst_markers,
        layout=layout
    )



def register_callbacks(app):
    @app.callback(
        Output('tree-graph-container', 'children'),
        [Input('upload-tree', 'contents'),
         Input('upload-metadata', 'contents'),
         Input('show-tip-labels', 'value')],
        [State('upload-tree', 'filename'),
         State('upload-metadata', 'filename')]
    )
    def update_tree_tab1(tree_contents, metadata_contents, show_labels, tree_filename, metadata_filename):
        print("Triggered update_tree_tab1 callback")  # Debug
        if tree_contents and metadata_contents:
            try:
                # Decode tree and metadata files
                tree_data = base64.b64decode(tree_contents.split(",", 1)[1]).decode("utf-8")
                metadata_data = base64.b64decode(metadata_contents.split(",", 1)[1]).decode("utf-8")
                tree_file = "uploaded_tree.tree"
                metadata_file = "uploaded_metadata.tsv"

                with open(tree_file, "w") as f:
                    f.write(tree_data)
                with open(metadata_file, "w") as f:
                    f.write(metadata_data)
                
                print(f"Tree and metadata files saved: {tree_filename}, {metadata_filename}")  # Debug

                show_tip_labels = 'SHOW' in show_labels

                # Generate tree plot
                fig = create_tree_plot(tree_file, metadata_file, show_tip_labels)
                print("Tree plot successfully created for Tab 1")  # Debug
                return dcc.Graph(figure=fig)
            except Exception as e:
                print(f"Error processing tree or metadata files in Tab 1: {str(e)}")  # Debug
                return html.Div(f"Error: {str(e)}", className="text-danger")
        print("Tree or metadata files missing for Tab 1")  # Debug
        return html.Div("Please upload both a tree file and a metadata file.", className="text-warning")

    @app.callback(
        Output('tree-graph-container-2', 'children'),
        [Input('upload-tree-2', 'contents'),
         Input('upload-metadata-2', 'contents'),
         Input('show-tip-labels-2', 'value')],
        [State('upload-tree-2', 'filename'),
         State('upload-metadata-2', 'filename')]
    )
    def update_tree_tab2(tree_contents, metadata_contents, show_labels, tree_filename, metadata_filename):
        print("Triggered update_tree_tab2 callback")  # Debug
        if not tree_contents or not metadata_contents:
            print("Tree or metadata files missing for Tab 2")  # Debug
            return html.Div("Please upload both a tree file and a metadata file.", className="text-warning")

        try:
            tree_data = base64.b64decode(tree_contents.split(",", 1)[1]).decode("utf-8")
            metadata_data = base64.b64decode(metadata_contents.split(",", 1)[1]).decode("utf-8")
            tree_file = "tree_tab2.tree"
            metadata_file = "metadata_tab2.tsv"

            with open(tree_file, "w") as f:
                f.write(tree_data)
            with open(metadata_file, "w") as f:
                f.write(metadata_data)
            
            print(f"Tree and metadata files saved for Tab 2: {tree_filename}, {metadata_filename}")  # Debug

            show_tip_labels = 'SHOW' in show_labels
            # Use a smaller size for the tree on this specific tab
            fig = create_tree_plot(tree_file, metadata_file, show_tip_labels, height=600, width=600)

            print("Tree plot successfully created for Tab 2")  # Debug
            return dcc.Graph(figure=fig)
        except Exception as e:
            print(f"Error processing tree or metadata files in Tab 2: {str(e)}")  # Debug
            return html.Div(f"Error: {str(e)}", className="text-danger")


    @app.callback(
        [Output('phylo-city-marker-status', 'children'),
        Output('phylo-marker-lat', 'value'),
        Output('phylo-marker-lon', 'value')],
        [Input('phylo-marker-city-btn', 'n_clicks')],
        [State('phylo-marker-city', 'value')]
    )
    def find_phylo_city_coordinates(n_clicks, city_name):
        """Find coordinates for the entered city name when clicking 'Find Location'."""
        if not n_clicks or not city_name:
            raise PreventUpdate

        lat, lon, error_msg = get_city_coordinates(city_name)
        if lat and lon:
            return f"✅ Found: {lat}, {lon}", lat, lon
        else:
            return f"⚠️ Error: {error_msg}", None, None


    # Callback for Folium Map Display
    @app.callback(
        Output('phylo-map-container', 'children'),
        [Input('upload-geojson', 'contents'),
        Input('map-city', 'value'),
        Input('map-lat', 'value'),
        Input('map-lon', 'value'),
        Input('map-zoom', 'value'),
        Input('phylo-add-marker-btn', 'n_clicks')],
        [State('phylo-marker-name', 'value'),
        State('phylo-marker-city', 'value'),
        State('phylo-marker-lat', 'value'),
        State('phylo-marker-lon', 'value')]
    )
    def update_phylo_folium_map(geojson_contents, city_name, latitude, longitude, zoom, n_clicks, marker_name, marker_city, marker_lat, marker_lon):
        global MARKERS
        
        # ✅ Clear the markers when the window is opened (i.e., first callback execution)
        if ctx.triggered_id is None:
            MARKERS.clear()

        # ✅ Convert city name to lat/lon if provided
        if marker_city:
            marker_lat, marker_lon, error_msg = get_city_coordinates(marker_city)
            if error_msg:
                return html.Div(f"⚠️ Error: {error_msg}", className="text-danger")

        # ✅ Ensure valid marker data before adding
        if ctx.triggered_id == "phylo-add-marker-btn" and marker_name:
            if marker_lat is None or marker_lon is None:
                return html.Div("⚠️ Error: Provide either a city name or latitude/longitude for the marker.", className="text-danger")

            MARKERS.append({"name": marker_name, "lat": float(marker_lat), "lon": float(marker_lon)})

        # ✅ Decode GeoJSON if uploaded
        geojson_data = None
        if geojson_contents:
            try:
                content_type, content_string = geojson_contents.split(',')
                decoded = base64.b64decode(content_string).decode('utf-8')
                geojson_data = json.loads(decoded)
            except Exception as e:
                return html.Div(f"⚠️ Error parsing GeoJSON: {str(e)}", className="text-danger")

        # ✅ Generate updated Folium map
        folium_map_html = phylo_map.generate_folium_map(geojson_data, latitude, longitude, zoom, MARKERS)

        return html.Iframe(
            srcDoc=folium_map_html,
            width="100%",
            height="600px",
            style={"border": "none"}
        )




    #standalone map tab
    @app.callback(
        [Output('standalone-city-marker-status', 'children'),
        Output('standalone-marker-lat', 'value'),
        Output('standalone-marker-lon', 'value')],
        [Input('standalone-marker-city-btn', 'n_clicks')],
        [State('standalone-marker-city', 'value')]
    )
    def find_city_coordinates(n_clicks, city_name):
        """Find coordinates for the entered city name when clicking 'Find Location'."""
        if not n_clicks or not city_name:
            raise PreventUpdate

        lat, lon, error_msg = get_city_coordinates(city_name)
        if lat and lon:
            return f"✅ Found: {lat}, {lon}", lat, lon
        else:
            return f"⚠️ Error: {error_msg}", None, None


    @app.callback(
        Output('standalone-map-container', 'children'),
        [Input('upload-standalone-geojson', 'contents'),
        Input('standalone-map-zoom', 'value'),
        Input('standalone-add-marker-btn', 'n_clicks')],
        [State('upload-standalone-geojson', 'filename'),
        State('standalone-marker-name', 'value'),
        State('standalone-marker-city', 'value'),
        State('standalone-marker-lat', 'value'),
        State('standalone-marker-lon', 'value')]
    )
    def update_standalone_map(geojson_contents, zoom, marker_clicks, filename, marker_name, marker_city, marker_lat, marker_lon):
        global STANDALONE_MARKERS

        latitude, longitude = 40.650002, -73.949997  # Default: New York
        geojson_data = None

        # ✅ Clear standalone markers on page load
        if ctx.triggered_id is None:
            STANDALONE_MARKERS.clear()

        # ✅ Decode GeoJSON if uploaded
        if geojson_contents:
            try:
                content_type, content_string = geojson_contents.split(',')
                decoded = base64.b64decode(content_string).decode('utf-8')
                geojson_data = json.loads(decoded)
            except Exception as e:
                return html.Div(f"⚠️ Error parsing GeoJSON: {str(e)}", className="text-danger")

        # ✅ Add a new marker if button clicked
        if ctx.triggered_id == "standalone-add-marker-btn" and marker_name:
            if marker_city:
                # Convert city name to coordinates
                marker_lat, marker_lon, error_msg = get_city_coordinates(marker_city)
                if error_msg:
                    return html.Div(f"⚠️ Error: {error_msg}", className="text-danger")

            if marker_lat is None or marker_lon is None:
                return html.Div("⚠️ Error: Provide either a city name or latitude/longitude for the marker.", className="text-danger")

            STANDALONE_MARKERS.append({"name": marker_name, "lat": float(marker_lat), "lon": float(marker_lon)})

        # ✅ Generate updated Folium map
        standalone_map_html = phylo_map.generate_folium_map(geojson_data, latitude, longitude, zoom, STANDALONE_MARKERS)

        return html.Iframe(
            srcDoc=standalone_map_html,
            width="100%",
            height="600px",
            style={"border": "none"}
        )

    @app.callback(
        Output("download-svg", "data"),
        Input("download-svg-btn", "n_clicks"),
        [Input('show-tip-labels', 'value')],  # Capture show_tip_labels
        prevent_initial_call=True
    )
    def export_svg(n_clicks, show_tip_labels):
        """Efficiently exports the phylogenetic tree as an SVG file."""
        
        # Generate the tree figure with respect to the user selection
        tree_fig = create_tree_plot("uploaded_tree.tree", "uploaded_metadata.tsv", show_tip_labels=show_tip_labels)

        # Save as SVG in memory
        svg_io = io.BytesIO()
        tree_fig.write_image(svg_io, format="svg", engine="kaleido")  # Explicit engine choice

        return dcc.send_bytes(svg_io.getvalue(), filename="phylogenetic_tree.svg")

