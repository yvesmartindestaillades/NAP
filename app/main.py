import pandas as pd
import sys
import base64

import matplotlib.pyplot as plt
import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from my_layout import layout
import io

pd.options.mode.chained_assignment = None  # default='warn'
path = '/Users/ymdt/src/dreem_nap/'
sys.path.append(path)
from dreem_nap.study import *
import dreem_nap
import plotly.express as px  # (version 4.7.0 or higher)
import plotly.graph_objects as go
from dash import Dash, html, Input, Output, ctx, dcc, State


def main(path_to_data, app):

    # -- Import and clean data (importing csv into pandas)
    # df = pd.read_csv("intro_bees.csv")
    global state, study

    study = Study(path_to_data= path_to_data, 
        min_cov_bases=0,
        filter_by='sample')

    state = {
        'plot_type': 'mutation_histogram',
        'min_cov_bases': 1000,
        **{k:study.df[k].iloc[0] for k in ['sample','construct','section','cluster']},
        'experimental_variable': 'DMS_conc_mM',
        }
    print(state)

    df = study.df
    df.reset_index(inplace=True)
    print(df[:5])

    # ------------------------------------------------------------------------------
    # App layout
    app.layout = layout(study)

    @app.callback(
        Output(component_id='min_cov_bases', component_property='max'),
        Input(component_id='upload-data', component_property='contents'),
        State('upload-data', 'filename'))

    def update_data(contents, filename):
        print('Update data')
        if contents is not None:
            content_type, content_string = contents.split(',')

            decoded = base64.b64decode(content_string)
            try:
                if 'csv' in filename:
                    # Assume that the user uploaded a CSV file
                    study.set_df(pd.read_csv(
                        io.StringIO(decoded.decode('utf-8'))))
                elif 'xls' in filename:
                    # Assume that the user uploaded an excel file
                    study.set_df(pd.read_excel(io.BytesIO(decoded)))
            except:
                pass
        return int(study.df.worst_cov_bases.max())



    @app.callback(
        Output(component_id='my_plot', component_property='figure'),
        Input(component_id='generate_button', component_property='n_clicks'))
    def update_plot(n_clicks):
        print('Update plot')
        fig = getattr(study,state['plot_type'])(**{k:v for k,v in state.items() if k not in ['plot_type']},use_iplot=False)['fig']
        print(state)
        return fig


    @app.callback(
        Output("save_html_text", "data"),
        Input("save_html", "n_clicks"),
        State('my_plot', 'figure'),
        prevent_initial_call=True,
    )
    def save_html(n_clicks, fig):
        return dcc.send_string(go.Figure(data=fig['data'], layout=fig['layout']).to_html(), filename="my_plot.html")



    # Refresh all the dropdowns
    @app.callback(Output('sample', 'options'),
                Output('construct', 'options'),
                Output('section', 'options'),
                Output('cluster', 'options'),
                Input('plot_type', 'value'),
                Input('min_cov_bases', 'value'),
                Input('experimental_variable', 'value'),
                Input('sample', 'value'),
                Input('construct', 'value'),
                Input('section', 'value'),
                Input('cluster', 'value'),
                Input('library_item_1', 'value'),
                Input('library_item_1_label', 'children'),
                Input('library_item_2', 'value'),
                Input('library_item_2_label', 'children'),
                Input('library_item_3', 'value'),
                Input('library_item_3_label', 'children'))

    def refresh_state_filter_dropdown(plot_type,min_cov_bases,experimental_variable,sample, construct, section, cluster, library_item_1, library_item_1_label, library_item_2, library_item_2_label, library_item_3, library_item_3_label):
        print('Refresh state filter dropdown')
        global state
        state = {}
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
                for f in ['sample','construct','section','cluster']:
                    if state[f] is not None:
                        break
                    filters.append(f)        
        df = study.get_df(**{k:state[k] for k in filters})
        out = []
        attributes = ['sample','construct','section','cluster']
        for attr in attributes:
            out.append([{"label": str(s), "value": s} for s in df[attr].unique()])
        if len(df) == 0:
            out = [[]*len(attributes)]
        print(f"{state=}")

        return out
    

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
        print(out)
        return out



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

