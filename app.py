from dash import Dash
import dash_bootstrap_components as dbc
from layout import app_layout  # ✅ Includes new tabs
import callbacks 
import logging

logging.basicConfig(level=logging.DEBUG)

dbc_css = "https://cdn.jsdelivr.net/gh/AnnMarieW/dash-bootstrap-templates/dbc.min.css"

# Initialize app
app = Dash(__name__, external_stylesheets=[dbc.themes.VAPOR, dbc_css])

app.layout = app_layout

# ✅ Register callbacks
callbacks.register_callbacks(app)

server = app.server  # For deployment

if __name__ == '__main__':
    app.run_server(debug=True, port=8050)
