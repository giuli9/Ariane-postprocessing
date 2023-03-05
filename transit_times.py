#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 11:29:32 2023

@author: giuli9
"""
import xarray as xr
import numpy as np
import weighted
from plotly.subplots import make_subplots
import plotly.io as pio
import pandas as pd
import plotly.express as px

ks1=xr.open_dataset('Lions/ariane_positions_quantitative.nc')
ks2=xr.open_dataset('Tyrrhenian/ariane_positions_quantitative.nc')
ks3=xr.open_dataset('Sicily/ariane_positions_quantitative.nc')

# converting time in days as numpy array
lmax=4749
lmin=1828
final_age1=((ks1.final_age/np.timedelta64(1, 's'))/86400).values
final_age2=((ks2.final_age/np.timedelta64(1, 's'))/86400).values
final_age3=((ks3.final_age/np.timedelta64(1, 's'))/86400).values

#mean transport in Sv
final_transp1=((ks1.final_transp.values)/10**6)/(1+lmax-lmin) 
final_transp2=((ks2.final_transp.values)/10**6)/(1+lmax-lmin) 
final_transp3=((ks3.final_transp.values)/10**6)/(1+lmax-lmin) 

# pandas Dataframe
data1 = {'time (day)':final_age1, 'transport (Sv)':final_transp1}
df1 = pd.DataFrame(data1)
total_transp1=sum(final_transp1)

data2 = {'time (day)':final_age2, 'transport (Sv)':final_transp2}
df2 = pd.DataFrame(data2)
total_transp2=sum(final_transp2)

data3 = {'time (day)':final_age3, 'transport (Sv)':final_transp3}
df3 = pd.DataFrame(data3)
total_transp3=sum(final_transp3)

x1=round(weighted.median(df1['time (day)'],df1['transport (Sv)']))
y1=round(weighted.quantile(df1['time (day)'],df1['transport (Sv)'],0.90))

x2=round(weighted.median(df2['time (day)'],df2['transport (Sv)']))
y2=round(weighted.quantile(df2['time (day)'],df2['transport (Sv)'],0.90))

x3=round(weighted.median(df3['time (day)'],df3['transport (Sv)']))
y3=round(weighted.quantile(df3['time (day)'],df3['transport (Sv)'],0.90))

# bin dimension is about a year
bins=round(28494/365)

#plotting three panels in one figure
figures = [px.histogram(df1,x='time (day)',y='transport (Sv)',
                 labels={'time (day)':'time (day)', 'transport (Sv)':'transport (Sv)'},
                 nbins=bins,opacity=0.8, color_discrete_sequence=['indianred'],width=920                 
                 ),
           px.histogram(df2,x='time (day)',y='transport (Sv)',
                            labels={'time (day)':'time (day)', 'transport (Sv)':'transport (Sv)'},
                            nbins=bins,opacity=0.8, color_discrete_sequence=['indianred'],width=920                  
                            ),
           px.histogram(df3,x='time (day)',y='transport (Sv)',
                            labels={'time (day)':'time (day)', 'transport (Sv)':'transport (Sv)'},
                            nbins=bins,opacity=0.8, color_discrete_sequence=['indianred'],width=920                  
                            )]
fig = make_subplots(rows=len(figures), cols=1, x_title='Time (years)',
                    y_title='Transport (Sv)', vertical_spacing=0.1, subplot_titles=(u'<b>Gulf of Lions \u2192 Gibraltar Strait</b>',u'<b>Northern Tyrrhenian \u2192 Gibraltar Strait</b>',u'<b>Sicily Strait \u2192 Gibraltar Strait</b>')) 

for i, figure in enumerate(figures):
    for trace in range(len(figure["data"])):
        fig.append_trace(figure["data"][trace], row=i+1, col=1)

fig.update_xaxes(range=[0, 28494], row=1, col=1)       
fig.update_xaxes(range=[0, 28494], row=2, col=1)
fig.update_xaxes(range=[0, 28494], row=3, col=1)
        
fig['layout']['xaxis']['tickmode'] = 'array'
fig['layout']['xaxis']['tickvals'] = [4749, 9498, 14247, 18996, 23745, 28494]
fig['layout']['xaxis']['ticktext'] = ['13', '26', '39', '52', '65', '78']     

fig['layout']['xaxis2']['tickmode'] = 'array'
fig['layout']['xaxis2']['tickvals'] = [ 4749, 9498, 14247, 18996, 23745, 28494]
fig['layout']['xaxis2']['ticktext'] = [ '13', '26', '39', '52', '65', '78']

fig['layout']['xaxis3']['tickmode'] = 'array'
fig['layout']['xaxis3']['tickvals'] = [ 4749, 9498, 14247, 18996, 23745, 28494]
fig['layout']['xaxis3']['ticktext'] = [ '13', '26', '39', '52', '65', '78']

fig['layout']['yaxis2']['tickmode'] = 'array'
fig['layout']['yaxis2']['tickvals'] = [ 0, 0.0002, 0.0004, 0.0006, 0.0008, 0.001]
fig['layout']['yaxis2']['ticktext'] = [ '0', '2e-4', '4e-4', '6e-4', '8e-4', '1e-3']

fig['layout']['yaxis']['tickmode'] = 'array'
fig['layout']['yaxis']['tickvals'] = [ 0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12]
fig['layout']['yaxis']['ticktext'] = [ '0', '2e-2', '4e-2', '6e-2', '8e-2', '0.1', '0.12']

fig['layout']['yaxis3']['tickmode'] = 'array'
fig['layout']['yaxis3']['tickvals'] = [ 0, 0.002, 0.004, 0.006, 0.008, 0.01, 0.012]
fig['layout']['yaxis3']['ticktext'] = [ '0', '2e-3', '4e-3', '6e-3', '8e-3', '1e-2', '1.2e-2']

fig['layout']['bargap']=0.2

fig.add_annotation(
    dict(
        x=15000, y=0.11, # annotation point
        xref='x1', 
        yref='y1',
        text="<b>T50%% = %4.0f days (~ 5 years)</b>"%x1,
        showarrow=False,
        font_size=26
        ))
fig.add_annotation(
    dict(
        x=15000, y=0.09, # annotation point
        xref="x1", 
        yref="y1",
        text="<b>T90%% = %4.0f days (~ 18 years)</b>"%y1,
        showarrow=False,
        font_size=26))
fig.add_annotation(
    dict(
        x=15000, y=0.0009, # annotation point
        xref='x2', 
        yref='y2',
        text="<b>T50%% = %4.0f days (~ 6 years)</b>"%x2,
        showarrow=False,
        font_size=26
        ))
fig.add_annotation(
    dict(
        x=15000, y=0.00075, # annotation point
        xref='x2', 
        yref='y2',
        text="<b>T90%% = %4.0f days (~ 22 years)</b>"%y2,
        showarrow=False,
        font_size=26))
fig.add_annotation(
    dict(
        x=15000, y=0.011, # annotation point
        xref='x3', 
        yref='y3',
        text="<b>T50%% = %4.0f days (~ 8 years)</b>"%x3,
        showarrow=False,
        font_size=26))
fig.add_annotation(
    dict(
        x=15000, y=0.009, # annotation point
        xref='x3', 
        yref='y3',
        text="<b>T90%% = %4.0f days (~ 28 years)</b>"%y3,
        showarrow=False,
        font_size=26)
    )
fig.update_traces(xbins=dict( # bins used for histogram
        start=0.0,
        end=28494,
        size=(28494/78)
    ))

# adding the vertical line indicating the median
fig.add_vline(x=x1, line_dash = 'dash', line_color = 'black',row=1,col=1)
fig.add_vline(x=x2,line_dash = 'dash', line_color = 'black',row=2,col=1)
fig.add_vline(x=x3, line_dash = 'dash', line_color = 'black',row=3,col=1)

# specifying the notations location manually
fig['layout']['annotations'][0]['xref']='x1'
fig['layout']['annotations'][0]['xanchor']='left'
fig['layout']['annotations'][0]['font']['size']=32

fig['layout']['annotations'][1]['xref']='x2'
fig['layout']['annotations'][1]['xanchor']='left'
fig['layout']['annotations'][1]['font']['size']=32

fig['layout']['annotations'][2]['xref']='x3'
fig['layout']['annotations'][2]['xanchor']='left'
fig['layout']['annotations'][2]['font']['size']=32

fig['layout']['annotations'][4]['xshift']=-85

fig['layout']['annotations'][4]['font']['size']=28
fig['layout']['annotations'][3]['font']['size']=28

fig.update_layout(autosize=False, margin_l=150,margin_r=70)

fig.update_xaxes(showgrid=True, tickfont=dict(size=25))
fig.update_yaxes(tickfont=dict(size=25))

fig.show()
pio.write_image(fig,'Transit_times.jpg',scale=10,width=1080, height=1080) 

