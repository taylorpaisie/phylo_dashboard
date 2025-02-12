import folium
import json

def generate_folium_map(geojson_data=None, latitude=40.650002, longitude=-73.949997, zoom=6, markers=[]):
    """Generates a Folium map with optional markers."""
    
    attr = ('&copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a> '
            'contributors, &copy; <a href="http://cartodb.com/attributions">CartoDB</a>')
    
    m = folium.Map(location=[latitude, longitude], zoom_start=zoom, tiles="OpenStreetMap", attr=attr)

    # ✅ Add GeoJSON layer if provided
    if geojson_data:
        folium.GeoJson(geojson_data, style_function=lambda x: {
            'color': '#58bbff', 'fillColor': '#58bbff', 'fillOpacity': 0.25
        }).add_to(m)

    # ✅ Add default marker
    # folium.Marker(
    #     location=[latitude, longitude],
    #     popup="Selected Location",
    #     icon=folium.Icon(color="blue", icon="info-sign"),
    # ).add_to(m)

        # ✅ Add user-defined markers
    # for marker in markers:
    #     folium.CircleMarker(
    #         location=[marker["lat"], marker["lon"]],
    #         popup=marker["name"],
    #         icon=folium.Icon(color="red", icon="flag"),
    #     ).add_to(m)

    radium=30
    folium.CircleMarker(
        location=[latitude, longitude],
        popup="Selected Location",
        color="black",
        weight=1,
        fill_opacity=0.6,
        opacity=1,
        fill_color="red",
    ).add_to(m)

    # ✅ Add user-defined markers
    for marker in markers:
        folium.CircleMarker(
            location=[marker["lat"], marker["lon"]],
            popup=marker["name"],
            color="black",
            weight=1,
            fill_opacity=0.6,
            opacity=1,
            fill_color="green",
        ).add_to(m)


    return m._repr_html_()

def generate_standalone_map(geojson_data=None, latitude=30, longitude=-80, zoom=6):
    """Generates a completely independent Folium map for the standalone viewer with zoom support."""
    
    m = folium.Map(location=[latitude, longitude], zoom_start=zoom, tiles="CartoDB positron")

    # ✅ Add GeoJSON layer if provided
    if geojson_data:
        folium.GeoJson(geojson_data, style_function=lambda x: {
            'color': 'lightgray', 'fillColor': 'lightgray', 'fillOpacity': 0.4
        }).add_to(m)

    # ✅ Add a different marker style
    radium=30
    folium.CircleMarker(
        location=[latitude, longitude],
        popup="Selected Location",
        color="black",
        weight=1,
        fill_opacity=0.6,
        opacity=1,
        fill_color="red",
    ).add_to(m)

    return m._repr_html_()



