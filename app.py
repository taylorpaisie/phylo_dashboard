import dash
from dash import html, dcc
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
from ete3 import Tree, TreeStyle
import base64
import os

# Function to render the phylogenetic tree
def render_ete_tree(newick_file, output_file="tree.png"):
    """
    Render the phylogenetic tree using ete3 and save it as an image.
    """
    # Load the tree
    tree = Tree(newick_file)

    # Define tree style
    ts = TreeStyle()
    ts.show_leaf_name = True  # Show leaf names
    ts.show_branch_length = True  # Show branch lengths
    ts.show_branch_support = True  # Show branch support

    # Render the tree to an image
    tree.render(output_file, tree_style=ts)

# Dash App
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# Path to your Newick file
newick_file = "parsnp.filtered_polymorphic_sites.fasta.tree"
output_image = "tree.png"

# Render the tree to an image
render_ete_tree(newick_file, output_image)

# Encode the image for display in Dash
encoded_image = None
if os.path.exists(output_image):
    with open(output_image, "rb") as img_file:
        encoded_image = base64.b64encode(img_file.read()).decode("utf-8")

# App Layout
app.layout = dbc.Container(
    [
        html.H1("Phylogenetic Tree Viewer"),
        html.Div("The phylogenetic tree below was rendered using ete3."),
        html.Img(
            id="tree-image",
            src=f"data:image/png;base64,{encoded_image}" if encoded_image else "",
            style={"width": "100%"},
        ),
    ],
    fluid=True,
)

if __name__ == "__main__":
    app.run_server(debug=True)
