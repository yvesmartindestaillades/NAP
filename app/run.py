import pandas as pd
import sys
import base64

import matplotlib.pyplot as plt
import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import my_layout
import io

pd.options.mode.chained_assignment = None  # default='warn'
path = '/Users/ymdt/src/dreem_nap/'
sys.path.append(path)
from dreem_nap.study import *
import dreem_nap
import plotly.express as px  # (version 4.7.0 or higher)
import plotly.graph_objects as go
from dash import Dash, html, Input, Output, ctx, dcc, State


def print_state(state):
    for k, v in state.items():
        print(str(k)+' '*(max([len(i) for i in state.keys()]) + 10 - len(k)), v)



def main(path_to_data, app):

    # -- Import and clean data (importing csv into pandas)
    # df = pd.read_csv("intro_bees.csv")
    global state, study, loaded_study

    loaded_study, study = Study(), Study()

    if path_to_data is None:
        study = Study(path_to_data= path_to_data, 
                             min_cov_bases=0,
                             filter_by='sample')
        print('Loaded df', study.df)

    state = {
        'plot_type': 'mutation_histogram',
        }

    # ------------------------------------------------------------------------------
    # App layout
    app.layout = my_layout.layout(study)

    @app.callback(
        Output('min_cov_bases','max'),
        Output('sample', 'value'),
        Output('construct', 'value'),
        Output('section', 'value'),
        Output('cluster', 'value'),
        Output('dataset_label', 'children'),
        Input('upload-data', 'contents'),
        State('filter_by','value'),
        State('upload-data', 'filename'),
        prevent_initial_call=True,
        )
    def update_data(contents, filter_by, filename):
        global study, loaded_study

        print('Update data')
        if contents is not None:
            content_type, content_string = contents.split(',')
            decoded = base64.b64decode(content_string)
            try:
                if 'csv' in filename:
                    # Assume that the user uploaded a CSV file
                    reading_fun = pd.read_csv
                elif 'xls' in filename:
                    # Assume that the user uploaded an excel file
                    reading_fun = pd.read_excel
                df = reading_fun(io.StringIO(decoded.decode('utf-8')))
                loaded_study, study = Study(), Study()
                loaded_study.set_df(df)
                study.set_df(df, filter_by=filter_by)
                print('Loaded df', study.df)
                return int(study.df.worst_cov_bases.max()), study.df['sample'].iloc[0], study.df['construct'].iloc[0], study.df['section'].iloc[0], study.df['cluster'].iloc[0], 'Dataset: {}'.format(filename)

            except:
                print('There was an error processing this file.')
                return 0



    @app.callback(
        Output(component_id='my_plot', component_property='figure'),
        Input(component_id='generate_button', component_property='n_clicks'),
        prevent_initial_call=True
        )
    def update_plot(n_clicks):
        print('Update plot')
        fig = getattr(study,state['plot_type'])(**{k:v for k,v in state.items() if k not in ['plot_type','experimental_variable']},use_iplot=False)['fig']
        return fig


    @app.callback(
        Output("save_html_text", "data"),
        Input("save_html", "n_clicks"),
        State('my_plot', 'figure'),
        prevent_initial_call=True
    )
    def save_html(n_clicks, fig):
        return dcc.send_string(go.Figure(data=fig['data'], layout=fig['layout']).to_html(), filename="my_plot.html")



    # Refresh all the dropdowns
    @app.callback(
                Output('experimental_variable', 'options'),
                Output('library_attributes', 'options'),
                Output('sample', 'options'),
                Output('construct', 'options'),
                Output('section', 'options'),
                Output('cluster', 'options'),
                Input('plot_type', 'value'),
                Input('min_cov_bases', 'value'),
                Input('filter_by', 'value'),
                Input('experimental_variable', 'value'),
                Input('library_attributes', 'value'),
                Input('sample', 'value'),
                Input('construct', 'value'),
                Input('section', 'value'),
                Input('cluster', 'value'),
                Input('library_item_1', 'value'),
                Input('library_item_1_label', 'children'),
                Input('library_item_2', 'value'),
                Input('library_item_2_label', 'children'),
                Input('library_item_3', 'value'),
                Input('library_item_3_label', 'children'),
                Input('base_type', 'value'),
                Input('base_index_menu', 'value'),
                Input('select_residues', 'value'),
                Input('unique_sequence', 'value'),
                prevent_initial_call=True
                )
    def refresh_state_filter_dropdown(plot_type,min_cov_bases,filter_by,experimental_variable,library_attributes,sample, construct, section, cluster, \
        library_item_1, library_item_1_label, library_item_2, library_item_2_label, library_item_3, library_item_3_label, base_type, \
        base_index_menu, select_residues, unique_sequence):

        print('Refresh state filter dropdown')
        attributes = ['experimental_variable','library_attributes','sample','construct','section','cluster']
        if study.df is None:
            return len(attributes)*[[]]
        global state
        state = {'base_index':compute_base_index(base_index_menu, select_residues, unique_sequence)}
        for k in refresh_state_filter_dropdown.__code__.co_varnames[:refresh_state_filter_dropdown.__code__.co_argcount]:
            if not k.startswith('library_item'):
                state[k] = eval(k)
        for i in range(1,4):
            try:
                state[locals()[eval(f'library_item_{i}_label')]] = locals()[eval(f'library_item_{i}')]
            except:
                pass

        filters = ['min_cov_bases']
        if ctx.triggered_id is not None:
            state[ctx.triggered_id] = ctx.triggered[0]['value']
            if min_cov_bases is not None:
                for f in attributes[1:]:
                    if state[f] is not None:
                        break
                    filters.append(f)     
        study_loc = study if filter_by == 'sample' else study.filter_by_study()  
        df = study_loc.get_df(**{k:state[k] for k in filters})
        out = []
        for attr in attributes:
            if attr == 'experimental_variable':
                out.append([{"label": str(s), "value": s} for s in my_layout.maybe_a_library(study_loc.df)])# if type(study.df[s].iloc[0]) in [int,float]])
            elif attr == 'library_attributes':
                out.append([{"label": str(s), "value": s} for s in my_layout.maybe_a_library(study_loc.df)])
            else:
                out.append([{"label": str(s), "value": s} for s in df[attr].unique()])
        print(f"\nUpdate: {ctx.triggered_id} = {ctx.triggered[0]['value']} \n")
        print_state(state)
        return out
    

    def compute_base_index(base_index_menu, select_residues, unique_sequence):
        if base_index_menu == 'unique_sequence_field':
            return unique_sequence
        elif base_index_menu == 'index_range_field':
            return select_residues
        else:
            return None

    @app.callback(
        Output('library_item_1', 'options'),
        Output('library_item_1', 'value'),
        Output('library_item_1_div', 'style'),
        Output('library_item_1_label', 'children'),
        Output('library_item_2', 'options'),
        Output('library_item_2', 'value'),
        Output('library_item_2_div', 'style'),
        Output('library_item_2_label', 'children'),
        Output('library_item_3', 'options'),
        Output('library_item_3', 'value'),
        Output('library_item_3_div', 'style'),
        Output('library_item_3_label', 'children'),
        Input('library_attributes', 'value'))
    def display_library_item(value):
        out = []
        for i in range(1,4):
            if i<=len(value):
                out.append([{'label': str(s), 'value': s} for s in study.df[value[i-1]].unique()])
                out.append(study.df[value[i-1]].iloc[0])
                out.append({'width': "25%",'display': 'inline-block'})
                out.append(value[i-1])
            else:
                out.append([])
                out.append(None)
                out.append({'display': 'none'})
                out.append('')
        return out

    @app.callback(
        Output('unique_sequence_field','style'),
        Output('index_range_field','style'),
        Output('select_residues_field','style'),
        Input('base_index_menu','value')
        )
    def display_base_index(base_index_menu):
        if base_index_menu == 'unique_sequence_field':
            return {'display': 'inline-block', 'width': '75%'}, {'display': 'none'}, {'display': 'none'}
        if base_index_menu == 'index_range_field':
            return {'display': 'none'}, {'display': 'inline-block', 'width': '25%'},  {'display': 'inline-block', 'width': '60%'}
        return {'display': 'none'}, {'display': 'none'}, {'display': 'none'}
 

    @app.callback(
        Output('select_residues','options'),
        Output('select_residues','value'),
        Input('index_range_start','value'),
        Input('index_range_end','value'))
    def display_select_residues(index_range_start, index_range_end):
        if index_range_start is None or index_range_end is None:
            return [], None
        return [{'label': str(i), 'value': i} for i in range(index_range_start, index_range_end+1)], [i for i in range(index_range_start, index_range_end+1)]


#@app.callback(
#    [Output(component_id='my_plot', component_property='figure')],
#    [Input(component_id='construct_dropdown', component_property='value')]
#)
#def update_construct(construct):
#    return update_graph(construct)


#app = run_standalone_app(layout, callbacks, header_colors, __file__)
#server = app.server

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    app = Dash(__name__)
    main(path_to_data='/Users/ymdt/src/dreem_nap/ex/Lauren/lau.csv', app=app)
    app.run_server(debug=True)

