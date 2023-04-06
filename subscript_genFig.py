import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default='svg'
import config 


class Figure:  

    """
    purpose:
        generate figures for fitting analysis 
        
    input: 
        experimental and/or model data from fitNR function
        
    output: 
        savefig method stores the figure in ../ouptut/~
    
    comment:
        to do - generalise layout, x and y axis objects to methods to reduce
        the volume of code 
    """
    

    def __init__(self,rows=1,cols=1,shared_xaxes=True,shared_yaxes=True,subplot_titles='',legend=dict(x=0.65,y=1,font=dict(size=16))):
        pass



    def plotStruct(self,Q,expNR,expNR_err,modelQ,modelNR,inputLabels,title,row=1,col=1,showlegend=False):

        # initialise figure 
        self.fig = go.Figure()
        
        nFiles = len(Q)
        for i in range(nFiles):
    
            # plot experimental data points on top of the model trace
            self.fig.add_trace(go.Scatter(
                            x=Q.get(i),
                            y=expNR.get(i),
                            error_y=dict(
                                    type='data', # value of error bar given in data coordinates
                                    array=expNR_err.get(i),
                                    visible=True),
                            showlegend=True,
                            name=inputLabels.get(i),
                            mode='markers',
                            marker_symbol='circle',
                            marker_color=config.col_light[i],
                            marker=dict(size=20),
                            )
                        )
    
            # add model trace 
            self.fig.add_trace(go.Scatter(
                            x=modelQ.get(i),
                            y=modelNR.get(i),
                            showlegend=False,
                            name=inputLabels.get(i),
                            mode='lines',
                            line=dict(width=5,
                                      color=config.col_dark[i]),
                            )
                        )




        self.fig.update_layout(font=dict(family='Latex',
                                 size=24,
                                 color='black'
                                 ),
                                margin=dict(r=20,t=20,b=10),
                                autosize=False,
                                width=700,
                                height=500,
                                plot_bgcolor='white',
                                showlegend = True,
                                legend=dict(
                                            x=0.4,
                                            y=1,
                                            bgcolor='rgba(0,0,0,0)',
                                            ),
                                legend_font_size=16,
                                )
    
    
        self.fig.update_xaxes(title_text=r'$\Huge{Q \ \rm{[\mathring{A}^{-1}]}}$', 
                          range=[0,0.35],
                          showline= True,
                          linecolor= 'black',
                          linewidth=3,
                          ticks= 'inside',
                          mirror='allticks',
                          tickwidth=3,
                          tickcolor='black',
                          tickfont=dict(family='Latex',
                                                size=24,
                                                color='black',
                                                ),
                          type='linear',
                          title_standoff=40,
                          )

    
           
        self.fig.update_yaxes(title_text=r'$\Huge{\rm{R}}$', 
                          showline= True,
                          linecolor= 'black',
                          linewidth=3,
                          ticks= 'inside',
                          minor_ticks="inside",
                          mirror='allticks',
                          tickwidth=3,
                          tickcolor='black',
                          tickfont=dict(family='Latex',
                                                size=24,
                                                color='black',
                                                ),
                          title_standoff=20,
                          type='log',
                          #tickformat = ".1SI", #".1r",
                          exponentformat ="power",
                          nticks=5,
                          range=[-6.2,0.2] # log range by exponent: 10^0=1, 10^5=100000
                          )




    def plotExcess(self,contrastList, global_objective):
        
        rho_d_list = []
        for i in range(len(contrastList)):
    
            # store each fitted parameter (thickness left to vary, SLD constant) in list 
            rho_d_list.append(global_objective.parameters.varying_parameters()[i].value)  
    
        print("\n\nScattering excess (rho*d) values:")
        print(*rho_d_list, sep ='\n')
    

        # initialise figure 
        self.fig = go.Figure()
        
        self.fig.add_trace(go.Scatter(
                        x=[i for i in range(len(rho_d_list))],
                        y=rho_d_list,
                        showlegend=False,
                        mode='lines',
                        line=dict(width=5,
                                  # color=,
                        )
                    ))
        
        self.fig.update_layout(font=dict(family='Latex',
                                size=24,
                                color='black'
                                ),
                               margin=dict(r=20,t=20,b=10),
                               autosize=False,
                               width=700,
                               height=500,
                               plot_bgcolor='white',
                               showlegend = True,
                               legend=dict(
                                           x=0.4,
                                           y=1,
                                           bgcolor='rgba(0,0,0,0)',
                                           ),
                               legend_font_size=16,
                               )
   
   
        self.fig.update_xaxes(title_text=r'$\Huge{Contrast}$', 
                         showline= True,
                         linecolor= 'black',
                         linewidth=3,
                         ticks= 'inside',
                         mirror='allticks',
                         tickwidth=3,
                         tickcolor='black',
                         tickfont=dict(family='Latex',
                                               size=24,
                                               color='black',
                                               ),
                         type='linear',
                         title_standoff=20,
                         )

   
        self.fig.update_yaxes(title_text=r'$\Huge{\rho*d}$', 
                         showline= True,
                         linecolor= 'black',
                         linewidth=3,
                         ticks= 'inside',
                         #minor_ticks="inside",
                         mirror='allticks',
                         tickwidth=3,
                         tickcolor='black',
                         tickfont=dict(family='Latex',
                                               size=24,
                                               color='black',
                                               ),
                         type='linear',
                         title_standoff=20,
                         )
        

    def save_figure(self,title):
    
        # PNG, sizes are in pixels: e.g. 700 by 500
        outputDir = '../output/fit-'+title+'.png'
        self.fig.write_image(outputDir) # width=700*300,height=500*300,scale=1 
        
        # SVG, sizes are inches * dpi 
        # outputDir = '../output/fit-'+title+'.svg'
        # fig.write_image(outputDir,width=1*720,height=0.714*720,scale=1)


    def show_figure(self):
        self.fig.show()


