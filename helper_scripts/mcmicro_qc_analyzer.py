#!/usr//bin/env python3
# -*- coding: utf-8 -*-

"""
================================================================================
MCMICRO Quantification QC Analyzer (Dynamic Version)
================================================================================
Version: 3.0

Purpose:
  This script is a generic tool that performs a comprehensive QC analysis on
  quantification data from an MCMICRO pipeline. It receives all configuration
  via command-line arguments, making it reusable for any run.

Output:
  A single, self-contained HTML file containing interactive plots.
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import os
import glob
import json
import argparse
import sys

def load_and_prepare_data(csv_dir):
    """Loads all CSVs from a directory, combines them, and identifies columns."""
    print(f"Loading all CSV files from {csv_dir}...")
    csv_files = glob.glob(os.path.join(csv_dir, '*.csv'))
    if not csv_files:
        print(f"Error: No CSV files found in directory {csv_dir}", file=sys.stderr)
        return None, None, None

    df_list = []
    for f in csv_files:
        try:
            df_single = pd.read_csv(f)
            core_id = os.path.basename(f).replace('.csv', '')
            df_single['CoreID'] = core_id
            df_list.append(df_single)
            print(f"  - Loaded {os.path.basename(f)} ({len(df_single)} cells)")
        except Exception as e:
            print(f"Could not read {f}. Error: {e}", file=sys.stderr)
            
    if not df_list:
        print("Error: Failed to load any valid dataframes.", file=sys.stderr)
        return None, None, None

    full_df = pd.concat(df_list, ignore_index=True)
    print(f"Successfully combined {len(csv_files)} files into a single dataset with {len(full_df)} total cells.")

    shape_cols = ['CellID', 'X_centroid', 'Y_centroid', 'Area', 'MajorAxisLength',
                  'MinorAxisLength', 'Eccentricity', 'Solidity', 'Extent', 'Orientation', 'CoreID']
    marker_cols = [col for col in full_df.columns if col not in shape_cols]
    
    print(f"Found {len(marker_cols)} marker columns and {len(shape_cols)} shape/ID columns.")
    return full_df, marker_cols, shape_cols

def normalize_data(df, marker_cols, cofactor):
    """Applies arcsinh transformation to all marker columns."""
    print(f"Normalizing marker data using arcsinh(x / {cofactor})...")
    df_norm = df.copy()
    for col in marker_cols:
        df_norm[col] = np.arcsinh(df_norm[col] / cofactor)
    return df_norm

def generate_heatmap_data(df, marker_cols, core_ids):
    """Pre-calculates correlation matrices and returns data ready for JSON."""
    print("Pre-calculating correlation heatmaps...")
    heatmap_data_dict = {}

    corr_matrix_all = df[marker_cols].corr()
    heatmap_data_dict['All Cores'] = {
        'z': corr_matrix_all.values.tolist(),
        'x': corr_matrix_all.columns.tolist(),
        'y': corr_matrix_all.index.tolist()
    }

    for core in core_ids:
        df_core = df[df['CoreID'] == core]
        corr_matrix_core = df_core[marker_cols].corr()
        heatmap_data_dict[core] = {
            'z': corr_matrix_core.values.tolist(),
            'x': corr_matrix_core.columns.tolist(),
            'y': corr_matrix_core.index.tolist()
        }
        
    return json.dumps(heatmap_data_dict)


def create_html_report(df_all_data, marker_cols, pixel_size, output_path, source_dir):
    """Generates a fully interactive HTML report with dynamic plots."""
    print("Generating dynamic HTML report...")

    df_all_data['Area_um2'] = df_all_data['Area'] * (pixel_size ** 2)
    core_ids = df_all_data['CoreID'].unique().tolist()
    data_json = df_all_data.to_json(orient='records')
    
    heatmap_data_json = generate_heatmap_data(df_all_data, marker_cols, core_ids)

    html_string = f"""
    <html>
    <head>
        <title>MCMICRO QC Report</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
            body {{ font-family: sans-serif; margin: 2em; background-color: #f4f6f6; }}
            h1, h2 {{ color: #2c3e50; border-bottom: 2px solid #bdc3c7; padding-bottom: 10px; }}
            .plot-container {{ display: flex; flex-wrap: wrap; justify-content: center; }}
            .plot {{ margin: 10px; background-color: white; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
            .controls {{ margin-bottom: 20px; padding: 15px; background-color: white; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        </style>
    </head>
    <body>
        <h1>MCMICRO Quantification QC Report</h1>
        <p><strong>Source Directory:</strong> {source_dir}</p>

        <div class="controls">
            <label for="coreSelector">Select TMA Core:</label>
            <select id="coreSelector" onchange="updateAllPlots()">
                <option value="All Cores">All Cores</option>
                {''.join([f"<option value='{core}'>{core}</option>" for core in core_ids])}
            </select>
        </div>

        <h2>Cell Shape Analysis</h2>
        <div class="plot-container">
            <div id="plot-area" class="plot"></div>
            <div id="plot-solidity" class="plot"></div>
            <div id="plot-area-solidity" class="plot"></div>
        </div>

        <h2>Marker Correlation</h2>
        <div class="plot-container">
            <div id="plot-heatmap" class="plot"></div>
        </div>
        
        <h2>Dynamic Co-expression Analysis</h2>
        <div class="controls" style="display: flex; align-items: center;">
            <div style="margin-right: 20px;">
                <label for="marker1">X-Axis Marker:</label><br>
                <select id="marker1" onchange="updateCoexpressionPlot()">{ ''.join([f"<option value='{col}'>{col}</option>" for col in marker_cols])}</select>
            </div>
            <div>
                <label for="marker2">Y-Axis Marker:</label><br>
                <select id="marker2" onchange="updateCoexpressionPlot()">{ ''.join([f"<option value='{col}'>{col}</option>" for col in marker_cols])}</select>
            </div>
        </div>
        <div class="plot-container">
            <div id="coexpression-plot" class="plot"></div>
        </div>

        <h2>Single Marker Intensity Distributions</h2>
        <div class="plot-container" id="intensity-plots-container">
        </div>

    <script>
        const fullData = {data_json};
        const markerCols = {json.dumps(marker_cols)};
        const heatmapData = {heatmap_data_json};

        function getFilteredData() {{
            const selectedCore = document.getElementById('coreSelector').value;
            if (selectedCore === 'All Cores') {{
                return fullData;
            }}
            return fullData.filter(row => row.CoreID === selectedCore);
        }}

        function drawShapePlots(data) {{
            const area_um2 = data.map(row => row.Area_um2);
            const solidity = data.map(row => row.Solidity);
            Plotly.newPlot('plot-area', [{{ x: area_um2, type: 'histogram', nbinsx: 100 }}], {{ title: 'Cell Area Distribution (um2)', xaxis: {{title: 'Area (um2)'}}, yaxis: {{title: 'Count'}} }});
            Plotly.newPlot('plot-solidity', [{{ x: solidity, type: 'histogram', nbinsx: 100 }}], {{ title: 'Cell Solidity Distribution', xaxis: {{title: 'Solidity'}}, yaxis: {{title: 'Count'}} }});
            Plotly.newPlot('plot-area-solidity', [{{ x: area_um2, y: solidity, mode: 'markers', type: 'scattergl', marker: {{ opacity: 0.5 }} }}], {{ title: 'Area vs. Solidity', xaxis: {{title: 'Area (um2)'}}, yaxis: {{title: 'Solidity'}} }});
        }}

        function drawHeatmap(selectedCore) {{
            const coreHeatmapData = heatmapData[selectedCore];
            const trace = {{
                z: coreHeatmapData.z, x: coreHeatmapData.x, y: coreHeatmapData.y,
                type: 'heatmap', colorscale: 'RdBu_r', zmin: -1, zmax: 1
            }};
            Plotly.newPlot('plot-heatmap', [trace], {{title: `Marker Correlation (${{selectedCore}})`}});
        }}

        function drawIntensityPlots(data) {{
            const container = document.getElementById('intensity-plots-container');
            container.innerHTML = '';
            markerCols.forEach(marker => {{
                const intensityData = data.map(row => row[marker]);
                const plotDiv = document.createElement('div');
                plotDiv.className = 'plot';
                container.appendChild(plotDiv);
                Plotly.newPlot(plotDiv, [{{ x: intensityData, type: 'histogram', nbinsx: 100 }}], {{ title: `Intensity: ${{marker}}`, xaxis: {{title: 'arcsinh(Intensity)'}}, yaxis: {{title: 'Count'}} }});
            }});
        }}
        
        function updateCoexpressionPlot() {{
            const data = getFilteredData();
            const marker1 = document.getElementById('marker1').value;
            const marker2 = document.getElementById('marker2').value;
            const x_data = data.map(row => row[marker1]);
            const y_data = data.map(row => row[marker2]);
            const trace_contour = {{ 
                x: x_data, y: y_data, type: 'histogram2dcontour', 
                colorscale: 'Hot', showscale: false, name: 'density',
                nbinsx: 150, nbinsy: 150 
            }};
            const trace_scatter = {{ 
                x: x_data, y: y_data, mode: 'markers', type: 'scattergl', 
                marker: {{ opacity: 0.1, color: 'yellow', size: 2 }}, 
                name: 'cells' 
            }};
            const layout = {{ 
                title: `Co-expression: ${{marker1}} vs. ${{marker2}}`, 
                xaxis: {{ title: `arcsinh(${{marker1}})` }}, 
                yaxis: {{ title: `arcsinh(${{marker2}})` }}, 
                showlegend: false,
                plot_bgcolor: '#111', paper_bgcolor: '#000',
                font: {{ color: 'white' }}
            }};
            Plotly.newPlot('coexpression-plot', [trace_contour, trace_scatter], layout);
        }}

        function updateAllPlots() {{
            const selectedCore = document.getElementById('coreSelector').value;
            const data = getFilteredData();
            drawShapePlots(data);
            drawIntensityPlots(data);
            drawHeatmap(selectedCore);
            updateCoexpressionPlot();
        }}
        
        window.onload = function() {{
            document.getElementById('marker2').selectedIndex = 1;
            updateAllPlots();
        }};

    </script>
    </body>
    </html>
    """
    
    with open(output_path, 'w') as f:
        f.write(html_string)
    print("Report generation complete.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate a QC report from MCMICRO quantification data.')
    parser.add_argument('--input-dir', required=True, help='Path to the quantification DIRECTORY containing CSV files.')
    parser.add_argument('--output-path', required=True, help='Path for the final HTML report.')
    parser.add_argument('--pixel-size', type=float, required=True, help='The size of one pixel in microns (µm).')
    parser.add_argument('--cofactor', type=int, default=5, help='The cofactor for the arcsinh transformation.')
    args = parser.parse_args()
    
    df, marker_cols, shape_cols = load_and_prepare_data(args.input_dir)
    
    if df is not None:
        df_norm = normalize_data(df, marker_cols, args.cofactor)
        create_html_report(df_norm, marker_cols, args.pixel_size, args.output_path, args.input_dir)
