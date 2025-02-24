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
from io import StringIO
from dash import Input, Output, State, dcc, html, dash_table
from dash.exceptions import PreventUpdate
import dash_bio as dashbio
from dash import ctx
import phylo_map
from Bio import Phylo, SeqIO
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

def create_tree_plot(tree_file, metadata_file, show_tip_labels):
    # Load tree and metadata
    tree = Phylo.read(tree_file, 'newick')
    metadata = pd.read_csv(metadata_file, sep='\t')

    # Validate required columns
    if 'taxa' not in metadata.columns or 'location' not in metadata.columns:
        raise ValueError("Metadata file must contain 'taxa' and 'location' columns.")

    metadata['location'] = metadata['location'].fillna('Unknown')

    # Map metadata to colors
    location_colors = {
        loc: qualitative.Plotly[i % len(qualitative.Plotly)]
        for i, loc in enumerate(metadata['location'].unique())
    }

    # Generate x and y coordinates
    x_coords = {}
    y_coords = {}
    max_y = 0

    def assign_coordinates(clade, x_start=0, y_start=0):
        nonlocal max_y
        branch_length = clade.branch_length if clade.branch_length else 0.001
        if clade.is_terminal():
            x_coords[clade] = x_start + branch_length
            y_coords[clade] = y_start
            max_y = max(max_y, y_start)
            return y_start + 1
        else:
            y_positions = []
            for child in clade.clades:
                y_start = assign_coordinates(child, x_start + branch_length, y_start)
                y_positions.append(y_start - 1)
            x_coords[clade] = x_start + branch_length
            y_coords[clade] = sum(y_positions) / len(y_positions)
            return y_start

    assign_coordinates(tree.root)

    # Create line shapes for branches
    line_shapes = []
    for clade in tree.find_clades(order='level'):
        x_start = x_coords[clade]
        y_start = y_coords[clade]
        if clade.clades:
            y_positions = [y_coords[child] for child in clade.clades]
            line_shapes.append(dict(
                type='line',
                x0=x_start, y0=min(y_positions),
                x1=x_start, y1=max(y_positions),
                line=dict(color='black', width=1)
            ))
            for child in clade.clades:
                x_end = x_coords[child]
                y_end = y_coords[child]
                line_shapes.append(dict(
                    type='line',
                    x0=x_start, y0=y_end,
                    x1=x_end, y1=y_end,
                    line=dict(color='black', width=1)
                ))

    # Create scatter points for tips
    tip_markers = []
    seen_locations = set()
    for clade in tree.get_terminals():
        x = x_coords[clade]
        y = y_coords[clade]
        meta_row = metadata[metadata['taxa'] == clade.name]
        location = meta_row['location'].iloc[0] if not meta_row.empty else 'Unknown'
        color = location_colors.get(location, 'gray')

        show_legend = location not in seen_locations
        if show_legend:
            seen_locations.add(location)

        tip_markers.append(go.Scatter(
            x=[x],
            y=[y],
            mode='markers+text' if show_tip_labels else 'markers',
            marker=dict(size=10, color=color, line=dict(width=1, color='black')),
            name=location if show_legend else None,
            text=f"<b>{clade.name}</b><br>Location: {location}" if show_tip_labels else "",
            textposition="middle right",
            hoverinfo='text',
            showlegend=show_legend
        ))

    layout = go.Layout(
        title='Phylogenetic Tree with Metadata',
        xaxis=dict(title='Evolutionary Distance', showgrid=True, zeroline=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-1, max_y + 1]),
        shapes=line_shapes,
        height=800,
        legend=dict(title="Locations", orientation="h", y=-0.2),
    )

    return go.Figure(data=tip_markers, layout=layout)

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
            fig = create_tree_plot(tree_file, metadata_file, show_tip_labels)
            print("Tree plot successfully created for Tab 2")  # Debug
            return dcc.Graph(figure=fig)
        except Exception as e:
            print(f"Error processing tree or metadata files in Tab 2: {str(e)}")  # Debug
            return html.Div(f"Error: {str(e)}", className="text-danger")


    # Callback for Folium Map Display
    @app.callback(
        Output('phylo-map-container', 'children'),
        [Input('upload-geojson', 'contents'),
        Input('map-city', 'value'),  # Accept city name input
        Input('map-lat', 'value'),
        Input('map-lon', 'value'),
        Input('map-zoom', 'value'),
        Input('add-marker-btn', 'n_clicks')],
        [State('marker-name', 'value'),
        State('marker-lat', 'value'),
        State('marker-lon', 'value')]
    )
    def update_folium_map(geojson_contents, city_name, latitude, longitude, zoom, n_clicks, marker_name, marker_lat, marker_lon):
        global MARKERS  # Preserve added markers

        # Convert city name to lat/lon if provided
        if city_name:
            geocoded_lat, geocoded_lon, error_msg = get_city_coordinates(city_name)
            if geocoded_lat and geocoded_lon:
                latitude, longitude = geocoded_lat, geocoded_lon
            else:
                return html.Div(f"⚠️ Error: {error_msg}", className="text-danger")

        # Decode GeoJSON if uploaded
        geojson_data = None
        if geojson_contents:
            content_type, content_string = geojson_contents.split(',')
            decoded = base64.b64decode(content_string).decode('utf-8')
            geojson_data = json.loads(decoded)

        # Add a new marker if button clicked
        if ctx.triggered_id == "add-marker-btn":
            if marker_name:
                # If lat/lon for the marker is provided, use them. Otherwise, use city name lookup.
                if marker_lat is not None and marker_lon is not None:
                    marker_lat, marker_lon = float(marker_lat), float(marker_lon)
                elif city_name:
                    marker_lat, marker_lon = latitude, longitude  # Use city geolocation
                else:
                    return html.Div("⚠️ Error: Provide either a city name or latitude/longitude for the marker.", className="text-danger")

                MARKERS.append({"name": marker_name, "lat": marker_lat, "lon": marker_lon})

        # Generate the updated Folium map
        folium_map_html = phylo_map.generate_folium_map(geojson_data, latitude, longitude, zoom, MARKERS)

        return html.Iframe(
            srcDoc=folium_map_html,
            width="100%",
            height="1000px",
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
        """Update Standalone Folium map based on GeoJSON upload, city search, and markers."""

        global STANDALONE_MARKERS  # ✅ Declare the global variable

        latitude, longitude = 40.650002, -73.949997  # Default: New York
        geojson_data = None

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

