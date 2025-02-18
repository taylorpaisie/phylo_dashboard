from dash import dcc, html
import dash_bootstrap_components as dbc

# About Tab Content
about_content = dbc.Container([
    html.H2("About This App"),
    html.P(
        "This application is designed to help visualize phylogenetic trees and geographic distributions. "
        "It allows users to upload tree structures, metadata, and geographic data to create interactive visualizations."
    ),
    html.P(
        "The Phylogenetic Tree Visualization tab provides a way to analyze evolutionary relationships, "
        "while the Map Visualization tab helps in locating specimens geographically."
    ),
    
    html.P([
        "This application was created and maintained by Taylor K. Paisie. Check out my ",
        html.A("website", href="https://taylorpaisie.github.io/", target="_blank", className="text-primary fw-bold"),
        " or my ",
        html.A("Github Repositories", href="https://github.com/taylorpaisie", target="_blank", className="text-primary fw-bold"),
        "."
    ]),
    html.P(["If you have any questions or issues about the Phylogenetic Tree and Map Visualization Dashboard, please submit an issue on the ",
        html.A("Phylogenetic Tree and Map Visualization Dashboard Github Repository", href="https://github.com/taylorpaisie/phylo_dashboard", target="_blank", className="text-primary fw-bold"),
        "."
    ])

], className="mt-4")

# How to Use Tab Content
how_to_use_content = dbc.Container([
    html.H2("How to Use This App"),
    html.P("Follow these steps to use the application effectively:"),
    html.Ol([
        html.Li("Upload a Newick tree file in the Phylogenetic Tree Visualization tab."),
        html.Li("Upload a corresponding metadata file for additional information."),
        html.Li("Optionally, upload a GeoJSON file to add geographic context."),
        html.Li("Use the Map Visualization tab to explore locations and add custom markers."),
        html.Li("Toggle the tip labels to show or hide tree node names."),
    ]),
    html.P("For more details, refer to the documentation or contact support.")
], className="mt-4")

# Exporting the layouts
about_tab = dcc.Tab(label='About', children=[about_content])
how_to_use_tab = dcc.Tab(label='How to Use', children=[how_to_use_content])
