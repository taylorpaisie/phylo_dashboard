from dash import Dash
import dash_bootstrap_components as dbc
from layout import app_layout  # ✅ Includes new tabs
import callbacks 
import logging

logging.basicConfig(level=logging.DEBUG)

# Initialize app
app = Dash(__name__, external_stylesheets=[dbc.themes.VAPOR])

app.layout = app_layout

# ✅ Register callbacks
callbacks.register_callbacks(app)

server = app.server  # For deployment

if __name__ == '__main__':
    app.run_server(debug=True, port=8051)
