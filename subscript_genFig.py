import plotly.graph_objects as go
import plotly.io as pio
pio.kaleido.scope.default_format = "png"
import config


class Figure:

    """
    purpose:
        generate figures for fitting analysis

    input:
        experimental and/or model data from fitNR function via sampleInfo

    output:
        savefig method stores the figure in ../ouptut/~

    comment:
        to do - generalise layout, x and y axis objects to methods to reduce
        the volume of code
    """


    def __init__(self,rows=1,cols=1,shared_xaxes=True,shared_yaxes=True,subplot_titles='',legend=dict(x=0.65,y=1,font=dict(size=16))):
        pass



    def plotStruct(self,title,sampleInfo,row=1,col=1,showlegend=False):

        # initialise figure
        self.fig = go.Figure()

        # sample number index for indexing against the colour lists in config
        idx_colour = 0


        for Q_exp, R_exp, R_err_exp, label in zip(\
                sampleInfo['Q_exp'],sampleInfo['R_exp'],sampleInfo['R_err_exp'],\
                    sampleInfo['label'],\
                ):

            # plot experimental data points on top of the model trace
            self.fig.add_trace(go.Scatter(
                            x=Q_exp,
                            y=R_exp,
                            error_y=dict(
                                    type='data', # value of error bar given in data coordinates
                                    array=R_err_exp,
                                    visible=True),
                            showlegend=True,
                            name=label,

                            mode='markers',
                            opacity=config.opacity,
                            marker=dict(
                                        color=config.col_dark[idx_colour],
                                        symbol=config.marker_symbol[idx_colour],
                                        size=config.marker_size,
                                        line=dict(
                                            color='black',
                                            width=0.5,
                                            ),
                                        )
                        )
                )
            idx_colour += 1


        if config.plotOnlyExp == False:

            idx_colour = 0
            for Q_model, R_model, label in zip(\
                    sampleInfo['Q_model'],sampleInfo['R_model'],sampleInfo['label'],\
                    ):

                # add model trace
                self.fig.add_trace(go.Scatter(
                                x=Q_model,
                                y=R_model,
                                showlegend=False,
                                name=label,
                                mode='lines',
                                line=dict(width=config.line_width,
                                          color=config.col_light[idx_colour]
                                          ),
                                )
                            )

                # update index to move to next colour
                idx_colour += 1



        self.fig.update_layout(font=dict(family='Latex',
                                 size=26,
                                 color='black'
                                 ),
                                autosize=False,
                                width=700,
                                height=600,
                                margin=dict(l=20,r=20,t=20,b=20),
                                plot_bgcolor='white',
                                showlegend = True,
                                legend=dict(
                                            x=config.legend_x,
                                            y=config.legend_y,
                                            bgcolor='rgba(0,0,0,0)',
                                            title=config.legend_title,
                                            ),
                                legend_font_size=config.legend_fs,
                                )


        self.fig.update_xaxes(title_text=r'$\huge{ \rm{Q \ \left[\mathring{A}^{-1}\right] }}$',
                          range=config.xrange,
                          #automargin=False,
                          showline= True,
                          linecolor= 'black',
                          linewidth=3,
                          ticks= 'inside',
                          mirror='allticks',
                          tickwidth=3,
                          tickcolor='black',
                          tickfont=dict(family='Latex',
                                                size=20,
                                                color='black',
                                                ),
                          type='linear',
                          title_standoff=40,
                          )



        self.fig.update_yaxes(title_text=r'$\huge{\rm{R}}$',
                          #automargin=False,
                          showline= True,
                          linecolor= 'black',
                          linewidth=3,
                          ticks= 'inside',
                          minor_ticks="inside",
                          mirror='allticks',
                          tickwidth=3,
                          tickcolor='black',
                          tickfont=dict(family='Latex',
                                                size=20,
                                                color='black',
                                                ),
                          title_standoff=10,
                          type='log',
                          #tickformat = ".1SI", #".1r",
                          exponentformat ="power",
                          nticks=5,
                          range=config.yrange
                          )




    def plotExcess(self, sampleInfo, global_objective):

        # store each fitted parameter (thickness left to vary, SLD constant) in list
        rho_d_list = []
        for i in range(len(sampleInfo)):
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

    def reorder(self):

        self.fig.data = (self.fig.data[2],self.fig.data[1],self.fig.data[3],\
                        self.fig.data[0],\
                        )


    def save_figure(self,title):

        # define file name via prefix
        if config.plotOnlyExp == True:
            prefix = 'raw-'
        else:
            prefix = 'fit-'

        # PNG, sizes are in pixels: e.g. 700 by 500
        outputDir = '../output/'+prefix+title+'.png'
        self.fig.write_image(outputDir,scale=5) #, width=700*300,height=500*300,scale=1)

        #outputDir = '../output/'+prefix+title+'.svg'
        #self.fig.write_image(outputDir)
        # SVG, sizes are inches * dpi
        # outputDir = '../output/fit-'+title+'.svg'
        # fig.write_image(outputDir,width=1*720,height=0.714*720,scale=1)


    def show_figure(self):
        self.fig.show(renderer='png')
