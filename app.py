from Bio import Phylo
from dash import dcc, html, Input, Output
import dash
import pandas as pd
import plotly.graph_objs as go

# Define a function to extract x and y coordinates and draw the tree
def get_x_coordinates(tree):
    xcoords = tree.depths(unit_branch_lengths=True)
    return xcoords

def get_y_coordinates(tree, dist=1.3):
    maxheight = tree.count_terminals()
    ycoords = {leaf: maxheight - i * dist for i, leaf in enumerate(reversed(tree.get_terminals()))}

    def calc_row(clade):
        for subclade in clade:
            if subclade not in ycoords:
                calc_row(subclade)
        ycoords[clade] = (ycoords[clade.clades[0]] + ycoords[clade.clades[-1]]) / 2

    if tree.root.clades:
        calc_row(tree.root)
    return ycoords

def draw_clade(clade, x_start, line_shapes, x_coords, y_coords):
    x_curr = x_coords[clade]
    y_curr = y_coords[clade]

    line_shapes.append(dict(type='line', x0=x_start, y0=y_curr, x1=x_curr, y1=y_curr, line=dict(color='black', width=1)))

    if clade.clades:
        y_top = y_coords[clade.clades[0]]
        y_bot = y_coords[clade.clades[-1]]
        line_shapes.append(dict(type='line', x0=x_curr, y0=y_bot, x1=x_curr, y1=y_top, line=dict(color='black', width=1)))

        for child in clade:
            draw_clade(child, x_curr, line_shapes, x_coords, y_coords)

# Function to create the tree plot
def create_tree_plot(tree_file, metadata_file):
    tree = Phylo.read(tree_file, 'newick')
    x_coords = get_x_coordinates(tree)
    y_coords = get_y_coordinates(tree)

    metadata = pd.read_csv(metadata_file, sep='\t')
    location_colors = {loc: f"hsl({i * 360 / len(metadata['location'].unique())}, 70%, 50%)" for i, loc in enumerate(metadata['location'].unique())}

    line_shapes = []
    draw_clade(tree.root, 0, line_shapes, x_coords, y_coords)

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
        title='Phylogenetic Tree with Colored Tips by Location',
        xaxis=dict(title='Evolutionary Distance', showgrid=True, zeroline=False, range=[0, max(x_coords.values()) * 1.1]),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[min(y_coords.values()) - 1, max(y_coords.values()) + 1]),
        shapes=line_shapes,
        height=800
    )
    return go.Figure(data=tip_markers, layout=layout)

# Dash app setup
app = dash.Dash(__name__, external_stylesheets=['https://stackpath.bootstrapcdn.com/bootstrap/4.5.2/css/bootstrap.min.css'])
app.title = "Phylogenetic Tree Viewer"

# Create the tree figure
tree_file = "parsnp.filtered_polymorphic_sites.fasta.tree"
metadata_file = "metadata.tsv"
tree_fig = create_tree_plot(tree_file, metadata_file)

# App layout
app.layout = html.Div([
    html.Div([
        html.H1("Phylogenetic Tree Viewer", style={"textAlign": "center", "marginBottom": "20px"}),
    ], className="jumbotron"),
    html.Div([
        dcc.Graph(figure=tree_fig, id="tree-graph", style={"height": "75vh"}),
    ], className="container-fluid"),
    html.Div([
        html.Label("Highlight Strain:"),
        dcc.Input(id="search-strain", type="text", placeholder="Enter strain name", style={"width": "100%", "marginBottom": "10px"}),
        html.Div(id="strain-highlight-message", style={"color": "red"})
    ], className="container"),
])

# Callback to highlight a specific strain
@app.callback(
    Output("strain-highlight-message", "children"),
    Input("search-strain", "value")
)
def highlight_strain(strain_name):
    if strain_name:
        if strain_name in pd.read_csv(metadata_file, sep='\t')['strain'].values:
            return f"Strain {strain_name} is highlighted on the tree."
        return f"Strain {strain_name} not found in the data."
    return ""

if __name__ == "__main__":
    app.run_server(debug=True, port=8051)
