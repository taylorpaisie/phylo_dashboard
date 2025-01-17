from Bio import Phylo
from dash import dcc, html, Input, Output, State
import dash
import pandas as pd
import plotly.graph_objs as go
import base64

# Define a function to calculate rectangular coordinates for the tree
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

# Function to create the rectangular tree plot
def create_tree_plot(tree_file, metadata_file):
    tree = Phylo.read(tree_file, 'newick')
    x_coords, y_coords = get_rectangular_coordinates(tree)

    metadata = pd.read_csv(metadata_file, sep='\t')
    location_colors = {loc: f"hsl({i * 360 / len(metadata['location'].unique())}, 70%, 50%)" for i, loc in enumerate(metadata['location'].unique())}

    line_shapes = []
    draw_clade_rectangular(tree.root, 0, line_shapes, x_coords, y_coords)

    # Create scatter points for tips
    tip_markers = []
    for clade, x in x_coords.items():
        if clade.is_terminal():
            y = y_coords[clade]
            if clade.name in metadata['strain'].values:
                meta_row = metadata[metadata['strain'] == clade.name].iloc[0]
                color = location_colors.get(meta_row['location'], 'gray')
                tip_markers.append(go.Scatter(
                    x=[x],
                    y=[y],
                    mode='markers',
                    marker=dict(size=10, color=color, line=dict(width=1, color='black')),
                    name=meta_row['location'],
                    text=f"{clade.name}<br>Location: {meta_row['location']}<br>Date: {meta_row['date']}",
                    hoverinfo='text'
                ))

    layout = go.Layout(
        title='Phylogenetic Tree (Rectangular Layout)',
        xaxis=dict(title='Evolutionary Distance', showgrid=True, zeroline=False, range=[0, max(x_coords.values()) * 1.1]),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[min(y_coords.values()) - 1, max(y_coords.values()) + 1]),
        shapes=line_shapes,
        height=800
    )
    return go.Figure(data=tip_markers, layout=layout)

# Dash app setup
app = dash.Dash(__name__, external_stylesheets=['https://stackpath.bootstrapcdn.com/bootstrap/4.5.2/css/bootstrap.min.css'])
app.title = "Phylogenetic Tree Viewer"

# App layout
app.layout = html.Div([
    html.Div([
        html.H1("Phylogenetic Tree Viewer", style={"textAlign": "center", "marginBottom": "20px"}),
    ], className="jumbotron"),
    html.Div([
        html.Label("Upload Tree File (.tree):"),
        dcc.Upload(
            id="upload-tree",
            children=html.Div(["Drag and Drop or ", html.A("Select File")]),
            style={"width": "100%", "height": "60px", "lineHeight": "60px", "borderWidth": "1px", "borderStyle": "dashed", "borderRadius": "5px", "textAlign": "center", "marginBottom": "10px"},
            multiple=False
        ),
        html.Label("Upload Metadata File (.tsv):"),
        dcc.Upload(
            id="upload-metadata",
            children=html.Div(["Drag and Drop or ", html.A("Select File")]),
            style={"width": "100%", "height": "60px", "lineHeight": "60px", "borderWidth": "1px", "borderStyle": "dashed", "borderRadius": "5px", "textAlign": "center", "marginBottom": "10px"},
            multiple=False
        ),
        html.Div(id="upload-message", style={"color": "red"})
    ], className="container"),
    html.Div([
        dcc.Graph(id="tree-graph", style={"height": "75vh"}),
    ], className="container-fluid"),
])

# Callback to handle file uploads
@app.callback(
    [Output("tree-graph", "figure"), Output("upload-message", "children")],
    [Input("upload-tree", "contents"), Input("upload-metadata", "contents")],
    [State("upload-tree", "filename"), State("upload-metadata", "filename")]
)
def update_tree(tree_contents, metadata_contents, tree_filename, metadata_filename):
    if tree_contents and metadata_contents:
        try:
            tree_data = base64.b64decode(tree_contents.split(",", 1)[1]).decode("utf-8")
            metadata_data = base64.b64decode(metadata_contents.split(",", 1)[1]).decode("utf-8")
            tree_file = "uploaded_tree.tree"
            metadata_file = "uploaded_metadata.tsv"

            with open(tree_file, "w") as f:
                f.write(tree_data)

            with open(metadata_file, "w") as f:
                f.write(metadata_data)

            fig = create_tree_plot(tree_file, metadata_file)
            return fig, "Tree and metadata successfully uploaded and visualized."
        except Exception as e:
            return go.Figure(), f"An error occurred: {str(e)}"
    return go.Figure(), "Please upload both a tree file and a metadata file."

if __name__ == "__main__":
    app.run_server(debug=True, port=8051)
