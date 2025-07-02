import os
import argparse
import glob


import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from PIL import Image

''' Define program parameters '''
parser = argparse.ArgumentParser(
                                prog = "PWM-scores Visualisation",
                                description = "Visualisation of PWM-scores for your protein sequence",
                                epilog = "For more information, see https://github.com/KazanovLab/ProteaseSpecificityModels"
                                )

parser.add_argument('path_to_file_with_PWM_scores', type=str, help="Enter the whole path to your file with PWM-scores (XLSX-file)")
parser.add_argument('-f', '--frame', type=str, default=None, help="Enter the frame for visualisation of the sewquence part, for example: 67-77")
args = parser.parse_args()

# Visualising using heatmap plot #
def create_heatmap(vis_data, init_file):
    file_name = init_file.split('.xlsx')[0]
    vis_data.sort_values('MEROPS_code', inplace=True)
    vis_data = vis_data.iloc[::-1, :]
    heatmap_rows = vis_data['MEROPS_code'].tolist()
    heatmap_cols = [c for c in vis_data.columns if '.' in c]
    heatmap_data = vis_data[heatmap_cols].values
    custom_data = np.stack((vis_data['MEROPS_name'], vis_data['Protease_class']), axis=1)
    print(custom_data)

    fig = go.Figure(data = go.Heatmap(
                            z = heatmap_data,
                            x = heatmap_cols,
                            y = heatmap_rows,
                            customdata = custom_data,
                            text=heatmap_data,
                            texttemplate="%{text:.2f}",
                            colorbar={"title":'''Cleavage score''', "tickfont_size":12},
                            hovertemplate='<b>Position:</b> %{x}<br>' + \
                                          '<b>MEROPS code:</b> %{y}<br>' + \
                                          '<b>MEROPS name:</b> %{customdata[0]}<br>' + \
                                          '<b>Protease class:</b> %{customdata[1]}<br>' + \
                                          '<b>Cleavage score:</b> %{z:.2f}' + \
                                          '<extra></extra>'
                            )
                    )
    fig.layout.xaxis.title = ""
    fig.layout.yaxis.title = ""
    title = dict(text=f"Cleavage score for {file_name}", font=dict(size=18), x=0.6),
    fig.update_layout(xaxis = dict(automargin=True, tickangle=-30),
                      yaxis = dict(automargin=True),
                      margin = dict(t=50, b=50, l=50, r=50)
                      )
    print('\nSaving the heatmap in interactive mode...')
    fig.write_html(file_name + ".plotly.html", auto_open=True)
    print('Saving the heatmap in interactive mode is well done!')

    ''' Data reading '''
print('#'*40 + f" {args.path_to_file_with_PWM_scores} " + '#'*40)
data = pd.read_excel(args.path_to_file_with_PWM_scores, dtype={'num_AA':str})

if args.frame is None:
    subset = data.copy()
else:
    frame = list(map(int, args.frame.split('-')))
    selected_columns = list(data.columns)[:12] + [i for i in list(data.columns) if ('.' in i) and (frame[0] <= int(i.split('.')[1]) < frame[1])]
    subset = data[selected_columns].copy()
    
create_heatmap(subset, args.path_to_file_with_PWM_scores)
