import sys
path = '/Users/ymdt/src/dreem_nap/'
sys.path.append(path)

import dreem_nap
from dreem_nap.study import *
from dash import Dash, dcc, html, Input, Output



def maybe_a_library(df):
    return [s for s in df.columns if (s not in ['sample', 'construct', 'section', 'cluster', 'index', 'num_reads', 'Unnamed: 0', 'num_aligned', 'start', 'end', 'experimental_variable']
            and df[s].dtype in [int, float, str])
            and not s.startswith('deltaG') and not s.startswith('structure') and not s.startswith('poisson')
            and not s.endswith('_bases')
            and not s.startswith('skips')
            and not s.startswith('ROI')]


def layout(study):
    return html.Div([

        html.H1("New Analysis Pipeline with Dash",
                style={'text-align': 'center'}),

        html.Div(className='row', children=[
            html.Div(children=[
                html.Label(['Dataset:'], style={
                           'font-weight': 'bold', "text-align": "center"}, id='dataset_label'),
                dcc.Upload(
                    id='upload-data',
                    children=html.Div([
                        'Drag and Drop or ',
                        html.A('Select Files')
                    ]),
                    style={
                        'width': '100%',
                        'height': '33px',
                        'lineHeight': '33px',
                        'borderWidth': '1px',
                        'borderStyle': 'dashed',
                        'borderRadius': '5px',
                        'textAlign': 'center',
                    },
                    # Allow multiple files to be uploaded
                    multiple=False
                ),
                html.Div(id='output-data-upload')
            ], style=dict(width='25%')),

            html.Div(children=[
                html.Label(['Filter by'], style={
                           'font-weight': 'bold', "text-align": "center"}),
                dcc.Dropdown(
                    id='filter_by',
                    options=[
                        {'label': 'Sample', 'value': 'sample'},
                        {'label': 'Study', 'value': 'study'},
                    ],
                    value='sample',
                    searchable=False,
                    clearable=False,
                )
            ], style=dict(width='25%')),


            html.Div(children=[
                html.Label(['Plot Type'], style={
                           'font-weight': 'bold', "text-align": "center"}),
                dcc.Dropdown(
                    id='plot_type',
                    options=[
                        {'label': 'Mutation Histogram',
                            'value': 'mutation_histogram'},
                        {'label': 'DeltaG per Sample',
                         'value': 'deltaG_per_sample'},
                        {'label': 'Experimental Variable across Samples',
                         'value': 'variable_exp_across_samples'},

                    ],
                    value='mutation_histogram',
                    searchable=False,
                    clearable=False,
                )
            ], style=dict(width='25%')),

            html.Div(children=[
                html.Label(['Minimal base coverage'], style={
                           'font-weight': 'bold', "text-align": "center"}),
                dcc.Input(
                    id="min_cov_bases",
                    type='number',
                    placeholder="input type number",
                    style={
                        'width': '100%',
                                'height': '31px',
                                'lineHeight': '31px',
                                'borderRadius': '5px',
                                'textAlign': 'center',
                    },
                    value=0,
                    max=999999999,
                    debounce=True,
                )
            ], style=dict(width='25%')),
        ], style=dict(display='flex')),


        html.Div(className='row', children=[
            html.Div(children=[
                html.Label(['Experimental Variable'], style={
                    'font-weight': 'bold', "text-align": "center"}),
                dcc.Dropdown(id="experimental_variable",
                             #       options=[{"label": str(s), "value": s} for s in maybe_a_library(study.df) if type(study.df[s].iloc[0]) in [int,float]],
                             multi=False,
                             value='DMS_conc_mM',
                             style={'width': "100%"},
                             ),
            ], style=dict(width='25%')),
            html.Div(children=[
                html.Label(['Libraries'], style={
                    'font-weight': 'bold', "text-align": "center"}),
                dcc.Dropdown(id="library_attributes",
                             #    options=[{"label": str(s), "value": s} for s in maybe_a_library(study.df)],
                             multi=True,
                             value=[],
                             style={'width': "100%"},
                             ),
            ], style=dict(width='25%')),
            html.Div(children=[
                html.Label(["Library item 1"], style={
                    'font-weight': 'bold', "text-align": "center"}, id="library_item_1_label"),
                dcc.Dropdown(id="library_item_1",
                             options=[],
                             multi=True,
                             value=[],
                             style={'width': "100%"},
                             ),
            ],
                id='library_item_1_div',
                style=dict(width='25%')),
            html.Div(children=[
                html.Label(["Library item 2"], style={
                    'font-weight': 'bold', "text-align": "center"}, id='library_item_2_label'),
                dcc.Dropdown(id="library_item_2",
                             options=[],
                             multi=True,
                             value=[],
                             style={'width': "100%"},
                             ),
            ],
                id='library_item_2_div',
                style=dict(width='25%')),
            html.Div(children=[
                html.Label(["Library item 3"], style={
                           'font-weight': 'bold', "text-align": "center"}, id='library_item_3_label'),
                dcc.Dropdown(id="library_item_3",
                             options=[],
                             multi=True,
                             value=[],
                             style={'width': "100%"},
                             ),
            ],
                id='library_item_3_div',
                style=dict(width='25%')),
        ], style=dict(display='flex')),


        html.Div(className='row', children=[
            html.Div(children=[
                html.Label(['Sample'], style={
                           'font-weight': 'bold', "text-align": "center"}),
                dcc.Dropdown(id="sample",
                             #    options=[{"label": str(s), "value": s} for s in study.df['sample'].unique()],
                             multi=True,
                             # value=study.df['sample'].unique()[0],
                             style={'width': "100%"},
                             ),
            ], style=dict(width='25%')),
            html.Div(children=[
                html.Label(['Construct'], style={
                           'font-weight': 'bold', "text-align": "center"}),
                dcc.Dropdown(id="construct",
                             #   options=[{"label": str(s), "value": s} for s in study.df['construct'].unique()],
                             multi=True,
                             # value=study.df['construct'].unique()[0],
                             style={'width': "100%"},
                             )
            ], style=dict(width='25%')),
            html.Div(children=[
                html.Label(['Section'], style={
                           'font-weight': 'bold', "text-align": "center"}),
                dcc.Dropdown(id="section",
                             #     options=[{"label": str(s), "value": s} for s in study.df['section'].unique()],
                             multi=True,
                             # value=study.df['section'].unique()[0],
                             style={'width': "100%"},
                             ),
            ], style=dict(width='25%')),
            html.Div(children=[
                html.Label(['Cluster'], style={
                           'font-weight': 'bold', "text-align": "center"}),
                dcc.Dropdown(id="cluster",
                             #      options=[{"label": str(s), "value": s} for s in study.df['cluster'].unique()],
                             multi=True,
                             # value=study.df['cluster'].unique()[0],
                             style={'width': "100%"},
                             ),
            ],
                style=dict(width='25%')),
        ], style=dict(display='flex')),



        html.Div(className='row', children=[
            html.Div(children=[
                html.Label(['Base type'], style={
                    'font-weight': 'bold', "text-align": "center"}),
                dcc.Checklist(
                    ['A', 'C', 'G', 'T'],
                    id='base_type',
                    value=['A', 'C', 'G', 'T'],
                    style={'width': "100%"},
                )
            ],
                style=dict(width='10%')),

            html.Div(children=[
                html.Label(['Index selection'], style={
                           'font-weight': 'bold', "text-align": "center"}),
                dcc.Dropdown(id="base_index_menu",
                             options=[{"label": "Unique sequence", "value": "unique_sequence_field"}, {
                                 'label': 'Range ', 'value': 'index_range_field'}],
                             multi=False,
                             style={'width': "100%"},
                             ),
            ],
                style={'width':'10%', 'margin-right':'10px'}),

            html.Div(children=[
                html.Label(['Unique sequence'], style={
                           'font-weight': 'bold', "text-align": "center"}),
                dcc.Input(id="unique_sequence",
                          type='text',
                          value='',
                          style={'width': "100%",'margin-right': "15px"},
                          debounce=True,
                          ),
            ],
                id='unique_sequence_field',
                style={'width': '25%'}),

            html.Div(children=[
                html.Label(['Index range'], style={
                           'font-weight': 'bold', "text-align": "center"}),
                html.Div(children=[
                    html.Label(['Start'], style={
                        'font-weight': 'bold', "text-align": "center"}),
                    dcc.Input(id="index_range_start",
                              type='number',
                              style={
                                  'width': "50px", 'margin-right': "15px", 'margin-left': "5px"},
                              ),
                    html.Label(['End'], style={
                        'font-weight': 'bold', "text-align": "center"}),
                    dcc.Input(id="index_range_end",
                              type='number',
                              style={
                                  'width': "50px", 'margin-right': "15px", 'margin-left': "5px"},
                              ),
                ], style=dict(display='flex')),
            ],
                id='index_range_field',
                style=dict(width='flex')),


            html.Div(children=[
                html.Label(['Select residues'], style={
                    'font-weight': 'bold', "text-align": "center"}),
                dcc.Dropdown(id="select_residues",
                             multi=True,
                             style={'width': "100%"},
                             ), ], 

                             style=dict(width='75%',display='none'),
                             id = 'select_residues_field'
                             ),


        ],
            style=dict(display='flex')),


        html.Div(className='row', children=[
            html.Div(children=[
                html.Button('Submit!', id='generate_button', n_clicks=0),
            ], 
                style=dict(width='flex')),
            html.Div(children=[
                html.Div([
                    html.Button('Save html', id='save_html', n_clicks=0),
                    dcc.Download(id="save_html_text")
                ])]),

        ], 
            style=dict(display='flex')),

        html.Br(),

        dcc.Graph(id='my_plot', figure={}),




        # Placeholder
        html.Button(id='placeholder1',  style={'display': 'none'}),
        html.Button(id='placeholder2',  style={'display': 'none'}),
        html.Button(id='placeholder3',  style={'display': 'none'}),
    ])
