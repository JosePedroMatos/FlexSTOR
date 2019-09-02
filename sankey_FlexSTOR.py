'''
Created on 28 Feb 2017

@author: Jose Pedro Matos

Produces Sankey plots of the KWO system.
'''

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

#===============================================================================
# from dateutil.relativedelta import relativedelta
#===============================================================================
from general import loadDailySystemData, loadEnergyPerM3
from matplotlib.sankey import Sankey as SK

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

def edge(func):
    def func_wrapper(*args, **kwargs):
        if 'edgecolor' not in kwargs.keys() and 'ec' not in kwargs.keys():
            kwargs['edgecolor'] = [0, 0, 0, 1]
        return func(*args, **kwargs)
    return func_wrapper

class SKC(SK):
    def __init__(self, details=True, *args, **kwargs):
        super(SKC, self).__init__(*args, **kwargs)
        
        self.details = details
        
    @edge
    def add(self, *args, **kwargs):
        for i0 in range(len(kwargs['flows'])):
            if abs(kwargs['flows'][i0])<=1E-6:
                if '-' in str(kwargs['flows'][i0])[0]:
                    kwargs['flows'][i0] = -1E-6
                else:
                    kwargs['flows'][i0] = 1E-6
                kwargs['labels'][i0] = None
        if self.details and 'label' in kwargs:
            if kwargs['patchlabel'] != None:
                kwargs['patchlabel'] = kwargs['label']
            else:
                kwargs['patchlabel'] = ''
        else:
            kwargs['labels'] = [None]*len(kwargs['labels'])
        super().add(*args, **kwargs)

def connection(sankey, inflow, outflow, where = 0.25, **kwargs):
    '''
    Creates secondary connections - a feature not available in the original matplotlib code.
    '''
    basePath = 0.25
    arrowEdge = sankey.shoulder
    tangent = 1/sankey.pitch
    tipExtra = arrowEdge/tangent
    strangeConstant0 = 0.08045018441136

    
    diagOut = sankey.diagrams[outflow[0]]
    diagIn = sankey.diagrams[inflow[0]]
    flowOut = diagOut.flows[outflow[1]]*sankey.scale
    flowIn = diagIn.flows[inflow[1]]*sankey.scale
    if abs(flowIn) + abs(flowOut)>sankey.scale*2E-6:
        xOut = diagOut.tips[outflow[1]][0]
        xIn = diagIn.tips[inflow[1]][0]
        yOut = diagOut.tips[outflow[1]][1]
        yIn = diagIn.tips[inflow[1]][1]
        angleOut = diagOut.angles[outflow[1]]%4
        angleIn = diagIn.angles[inflow[1]]%4
    
        # check distance
        x = xOut-xIn
        y = yOut-yIn
        
        if abs(x)<1E-6:
            # Perfectly aligned horizontally
            trunklength = abs(y)-tipExtra
            orientations = [0, 0]
            pathlengths = [basePath, basePath]
        elif abs(y)<1E-6:
            # Perfectly aligned vertically
            trunklength = abs(x)-tipExtra
            orientations = [0, 0]
            pathlengths = [basePath, basePath]
        elif (angleOut+angleIn)%2==0:
            # Same direction
            if angleOut%2==0:
                # Horizontal Out
                trunklength = abs(y)-(flowOut-flowIn)*0.5
                pathlength = abs(x)-tipExtra+flowIn*strangeConstant0-flowOut*(1-strangeConstant0)
                pathlengths = [pathlength*where, pathlength*(1-where)]
                if y*x<0:
                    orientations = [-1, 1]
                else:
                    orientations = [1, -1]    
                pass
            else:
                # Vertical Out
                trunklength = abs(x)-(flowOut-flowIn)*0.5
                pathlength = abs(y)-tipExtra+flowIn*strangeConstant0-flowOut*(1-strangeConstant0)
                pathlengths = [pathlength*where, pathlength*(1-where)]
                if y*x<0:
                    orientations = [1, -1]
                else:
                    orientations = [-1, 1]        
        else:
            # Other direction
            if angleOut%2==0:
                # Horizontal Out
                trunklength = abs(x)-flowOut*(-strangeConstant0+0.5)-tipExtra+flowIn*0.5
                pathlengths = [abs(y)+flowIn*strangeConstant0, basePath]
                if x*y>0: 
                    orientations = [-1, 0]
                else:
                    orientations = [1, 0]
            else:
                # Vertical Out
                trunklength = abs(y)-flowOut*(-strangeConstant0+0.5)-tipExtra+flowIn*0.5
                pathlengths = [abs(x)+flowIn*strangeConstant0, basePath]
                if x*y>0: 
                    orientations = [1, 0]
                else:
                    orientations = [-1, 0]
            
        flows = [-flowIn/sankey.scale, -flowOut/sankey.scale]
        labels = [None, None]
        
        if 'facecolor' not in kwargs.keys():
            kwargs['facecolor'] = diagIn.patch.get_facecolor()
        if 'edgecolor' not in kwargs.keys():
            kwargs['edgecolor'] = diagIn.patch.get_edgecolor()
        if 'patchlabel' not in kwargs.keys():
            kwargs['patchlabel'] = ''
        
        sankey.add(flows=flows, orientations=orientations, labels=labels, prior=outflow[0], connect=(outflow[1], 1), trunklength=trunklength, pathlengths=pathlengths, **kwargs)

def KWOSankey(a, block=True, title='', scale=1, units='Mm³/day', color=None, palette='deep', ax=None, legend=True, grid=True, details=False, gridStep=100):
    '''
    Sankey plot customized for the KWO system
    '''
    
    matplotlib.rcParams.update({'font.size': 10, 'font.weight': 'bold', 'patch.linewidth': 0.6, 'legend.fontsize': 10})
    
    names = ['Aare', 'Innertkirchen 1', 'Handeck', 'Innertkirchen 2', 'Hopflauenen',
             'Leimboden', 'Trift', 'Fuhren', 'Handeck 2', 'Handeck 1', 
             'Räterichsboden', 'Handeck 3', 'Gelmer', 'Grimsel', 'Oberaar', 
             '']    
    if color==None:
        palette = sns.color_palette(palette, n_colors=len(names))
    else:
        palette = sns.light_palette(color, n_colors=len(names), reverse=True)
    sns.set_palette(palette, n_colors=len(names))
    colors = iter(palette)
    next(colors)
    
    if ax==None:
        fig = plt.figure(figsize=cm2inch(15, 13), dpi=80)
        ax = fig.add_subplot(1, 1, 1)
    
    # Transformation to Mm3/day (Important)
    a=a*86400/1E6
    
    sankey = SKC(ax=ax, format='%.3G', unit=units, scale=scale, head_angle=100, shoulder=0.03, offset=0.2, details=details)

    fuhrenToTrift = a['Pum_Fuhren [m3/s]'] + a['Tur_Fuhren [m3/s]']
    triftToInnertkirchen = a['Zufl_Trift [m3/s]'] + fuhrenToTrift + a['Tur_Han3IsoT [m3/s]'] + a['Pum_Han3Dipu [m3/s]'] - a['Tur_HopfTrift [m3/s]'] - a['Pum_Han3IsoT [m3/s]']
    if triftToInnertkirchen<0:
        triftToInnertkirchen = 0
    handeckInMinusTrift = a['Zufl_Handeck [m3/s]']-triftToInnertkirchen
    handeckToInnertkirchen = a['Tur_Han1 [m3/s]'] + a['Tur_Han2 [m3/s]'] + a['Tur_Han3IsoH [m3/s]'] + handeckInMinusTrift - a['Pum_Han3Dipu [m3/s]']-a['Pum_Han3IsoH [m3/s]']
    hopflauenenToInnertkirchen = a['Tur_HopfTrift [m3/s]'] + a['Tur_HopfLeimb [m3/s]'] + a['Zufl_Hopf [m3/s]']
    grimselStorage = -(a['Zufl_Grimsel [m3/s]']+a['Tur_Gri2 [m3/s]']-a['Tur_Gri1M2 [m3/s]']-a['Pum_Gri2 [m3/s]']-a['Zuleit_Gri-Gel [m3/s]'])
    oberaarStorage = -(a['Zufl_Oberaar [m3/s]']+a['Pum_Gri2 [m3/s]']-a['Tur_Gri2 [m3/s]']-a['Tur_Gri1M1 [m3/s]'])
    raboStorage = -(a['Zufl_Räbo [m3/s]']+a['Tur_Gri1M2 [m3/s]']+a['Tur_Gri1M1 [m3/s]']-a['Tur_Han2 [m3/s]']-a['Tur_Han3IsoT [m3/s]']-a['Tur_Han3IsoH [m3/s]']+a['Pum_Han3IsoT [m3/s]']+a['Pum_Han3IsoH [m3/s]'])
    gelmerStorage = -(a['Zufl_Gelmer [m3/s]']+a['Zuleit_Gri-Gel [m3/s]']-a['Tur_Han1 [m3/s]'])

    # 0.Aare
    flows = [a['Tur_Inn1 [m3/s]'], a['Tur_Inn2 [m3/s]'], -a['Tur_Inn1 [m3/s]']-a['Tur_Inn2 [m3/s]']]
    labels = [None, None, ' ']
    orientations = [1, -1, 0]
    pathlengths = [1.5, 1.5, 0.25]
    sankey.details = True
    sankey.add(flows=flows, orientations=orientations, labels=labels, rotation=-90, patchlabel=None, pathlengths=pathlengths, label='Aare', facecolor=next(colors))
    sankey.details = details
    
    # 1.Innertkirchen 1
    flows = [triftToInnertkirchen, handeckToInnertkirchen, -a['Tur_Inn1 [m3/s]']]
    labels = ['', '', None]
    orientations = [-1, 0, -1]
    pathlengths = [0.25, 0.25, 1]
    sankey.add(flows=flows, orientations=orientations, labels=labels, prior=0, connect=(0, 2), patchlabel='', trunklength=1, pathlengths=pathlengths, label='Innertkirchen 1', facecolor=next(colors))
    
    # 2.Handeck
    flows = [handeckInMinusTrift, a['Tur_Han1 [m3/s]'], a['Tur_Han2 [m3/s]'], a['Tur_Han3IsoH [m3/s]'], -a['Pum_Han3Dipu [m3/s]'], -handeckToInnertkirchen, -a['Pum_Han3IsoH [m3/s]']]
    labels = ['', '', '', '', '', None, '']
    orientations = [1, 1, 0, -1, -1, 0, -1]
    pathlengths = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]
    sankey.add(flows=flows, orientations=orientations, labels=labels, prior=1, connect=(1, 5), patchlabel='', trunklength=1, pathlengths=pathlengths, label='Handeck', facecolor=next(colors))
    
    # 3.Innertkirchen 2
    flows = [a['Zufl_Hopf [m3/s]'], hopflauenenToInnertkirchen-a['Zufl_Hopf [m3/s]'], -a['Tur_Inn2 [m3/s]']]
    labels = ['', '', '']
    orientations = [-1, 0, 1]
    pathlengths = [0.75, 0.25, 1.5]
    sankey.add(flows=flows, orientations=orientations, labels=labels, prior=0, connect=(1, 2), patchlabel='', trunklength=0.5, pathlengths=pathlengths, label='Innertkirchen 2', facecolor=next(colors))
    
    # 4.Hopflauenen
    flows = [a['Tur_HopfTrift [m3/s]'], a['Tur_HopfLeimb [m3/s]'], -hopflauenenToInnertkirchen + a['Zufl_Hopf [m3/s]']]
    labels = ['', None, None]
    orientations = [0, -1, 0]
    pathlengths = [0.25, 0.25, 0.25]
    sankey.add(flows=flows, orientations=orientations, labels=labels, prior=3, connect=(1, 2), patchlabel='', trunklength=1.5, pathlengths=pathlengths, label='Hopflauenen', facecolor=next(colors))
    
    # 5.Leimboden
    flows = [a['Zufl_Leimb [m3/s]'], -a['Tur_HopfLeimb [m3/s]']]
    labels = ['', None]
    orientations = [1, 0]
    sankey.add(flows=flows, orientations=orientations, labels=labels, prior=4, connect=(1, 1), patchlabel='', trunklength=1.5, label='Leimboden', facecolor=next(colors))
    
    # 6.Trift
    flows = [a['Pum_Han3Dipu [m3/s]'], a['Tur_Han3IsoT [m3/s]'], a['Zufl_Trift [m3/s]'], fuhrenToTrift, -a['Tur_HopfTrift [m3/s]'], -a['Pum_Han3IsoT [m3/s]'], -triftToInnertkirchen]
    labels = ['', '','', None, None, '', None]
    orientations = [1, 1, 0, -1, 0, 1,  1]
    pathlengths = [1, 1.5, 0.25, 0.25, 1.5, 0.25, 0.25]
    sankey.add(flows=flows, orientations=orientations, labels=labels, prior=4, connect=(0, 4), patchlabel='', trunklength=0.5, pathlengths=pathlengths, label='Trift', facecolor=next(colors))
    
    # 7.Fuhren
    flows = [a['Zufl_Teuflaui [m3/s]'], a['Zufl_Fuhren [m3/s]'], -fuhrenToTrift]
    labels = ['', '' , None]
    orientations = [1, -1, 0]
    pathlengths=[0.5, 0.25, 1]
    sankey.add(flows=flows, orientations=orientations, labels=labels, prior=6, connect=(3, 2), patchlabel='', trunklength=0.7, pathlengths=pathlengths, label='Fuhren', facecolor=next(colors))
        
    # 8.Handeck 2
    flows = [a['Tur_Han2 [m3/s]'], -a['Tur_Han2 [m3/s]']]
    labels = [None]*len(flows)
    orientations = [0, 0]
    pathlengths=[0.25, 0.25]
    sankey.add(flows=flows, orientations=orientations, labels=labels, prior=2, connect=(2, 1), patchlabel='', trunklength=0.75, pathlengths=pathlengths, label='Handeck 2', facecolor=next(colors))
    
    # 9.Handeck 1
    flows = [a['Tur_Han1 [m3/s]'], -a['Tur_Han1 [m3/s]']]
    labels = ['' , None]
    orientations = [0, -1]
    pathlengths=[0.25, 7]
    sankey.add(flows=flows, orientations=orientations, labels=labels, prior=2, connect=(1, 1), patchlabel='', trunklength=0.75, pathlengths=pathlengths, label='Handeck 1', facecolor=next(colors))
    
    # 10. Räbo
    flows = [a['Zufl_Räbo [m3/s]'], a['Tur_Gri1M2 [m3/s]'], a['Tur_Gri1M1 [m3/s]'], -a['Tur_Han2 [m3/s]'], -a['Tur_Han3IsoT [m3/s]'], -a['Tur_Han3IsoH [m3/s]'], a['Pum_Han3IsoT [m3/s]'], a['Pum_Han3IsoH [m3/s]'], min(-1E-7, raboStorage), max(1E-7, raboStorage)]
    labels = ['', '', '' , None, '' , '', '', '', 'Filling', 'Emptying']
    orientations = [1, 1, 0, 0, -1, -1, -1, -1, 1, 1]
    pathlengths=[0.25, 0.5, 0.5, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]
    sankey.add(flows=flows, orientations=orientations, labels=labels, prior=8, connect=(0, 3), patchlabel='', trunklength=1, pathlengths=pathlengths, label='Räterichsboden', facecolor=next(colors))
    
    # 11.Handeck 3
    flows = [-a['Pum_Han3Dipu [m3/s]'], -a['Tur_Han3IsoT [m3/s]'], -a['Pum_Han3IsoT [m3/s]'], -a['Pum_Han3IsoH [m3/s]'], -a['Tur_Han3IsoH [m3/s]'], a['Pum_Han3IsoT [m3/s]'], a['Pum_Han3Dipu [m3/s]'], a['Pum_Han3IsoH [m3/s]'], a['Tur_Han3IsoT [m3/s]'], a['Tur_Han3IsoH [m3/s]']]
    labels = [None]*len(flows)
    orientations = [-1, -1, 1, 1, 1, -1, 1, 1, 1, 1]
    pathlengths=[0.25, 2.5, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]
    sankey.add(flows=flows, orientations=orientations, labels=labels, prior=6, connect=(1, 1), patchlabel='', trunklength=1.5, pathlengths=pathlengths, label='Handeck 3', facecolor=next(colors))

    # 12. Gelmer
    flows = [a['Zufl_Gelmer [m3/s]'], a['Zuleit_Gri-Gel [m3/s]'], -a['Tur_Han1 [m3/s]'],  min(-1E-7, gelmerStorage), max(1E-7, gelmerStorage)]
    labels = ['', '', None, 'Filling', 'Emptying']
    orientations = [1, 0, 0, 1, 1]
    pathlengths=[0.25, 0.25, 0.25, 0.25, 0.25]
    sankey.add(flows=flows, orientations=orientations, labels=labels, prior=9, connect=(0, 2), patchlabel='', trunklength=2, pathlengths=pathlengths, label='Gelmer', facecolor=next(colors))

    # 13. Grimsel
    flows = [a['Zufl_Grimsel [m3/s]'], -a['Tur_Gri1M2 [m3/s]'], a['Tur_Gri2 [m3/s]'], -a['Pum_Gri2 [m3/s]'], -a['Zuleit_Gri-Gel [m3/s]'], min(-1E-7, grimselStorage), max(1E-7, grimselStorage)]
    labels = ['', '', '', '',  None, 'Filling', 'Emptying']
    orientations = [1, -1, -1, -1, 0, 1, 1]
    pathlengths=[0.25, 0.25, 0.25, 0.25, 2, 0.25, 0.25]
    sankey.add(flows=flows, orientations=orientations, labels=labels, prior=12, connect=(1, 4), patchlabel='', trunklength=1.5, pathlengths=pathlengths, label='Grimsel', facecolor=next(colors))

    # 14. Oberaar
    flows = [a['Zufl_Oberaar [m3/s]'], a['Pum_Gri2 [m3/s]'], -a['Tur_Gri2 [m3/s]'], -a['Tur_Gri1M1 [m3/s]'], min(-1E-7, oberaarStorage), max(1E-7, oberaarStorage)]
    labels = ['', None, None, None, 'Filling', 'Emptying']
    orientations = [1, -1, -1, 1, 1, 1]
    pathlengths=[0.25, 0.25, 0.25, 0.25, 0.25, 0.25]
    sankey.add(flows=flows, orientations=orientations, labels=labels, prior=13, connect=(3, 1), patchlabel='', trunklength=1.5, pathlengths=pathlengths, label='Oberaar', facecolor=next(colors))

    # Trift Innertkirchen
    connection(sankey, inflow=(6, 6), outflow=(1, 0), where=0.65)
    
    # Grimsel Räbo
    connection(sankey, inflow=(13, 1), outflow=(10, 1), where=0.05)
    
    # Oberaar Räbo
    connection(sankey, inflow=(14, 3), outflow=(10, 2), where=0.5)
    
    # pIsoT in
    connection(sankey, inflow=(6, 5), outflow=(11, 5), where=0.05)#, patchlabel='pump IsoT')
     
    # Dipu out
    connection(sankey, inflow=(11, 0), outflow=(6, 0), where=0.5)#, patchlabel='Dipu')
     
    # IsoH out
    connection(sankey,  inflow=(11, 4), outflow=(2, 3), where=0.55)#, patchlabel='IsoH')
     
    # pIsoH out
    connection(sankey,  inflow=(11, 3), outflow=(10, 7), where=0.15)#, patchlabel='pumpIsoH')
     
    # IsoT out
    connection(sankey,  inflow=(11, 2), outflow=(10, 6), where=0.05)#, patchlabel='IsoT')
     
    # pIsoH in
    connection(sankey,  inflow=(2, 6), outflow=(11, 7), where=0.55)#, patchlabel='pump IsoH')
     
    # IsoT in
    connection(sankey,  inflow=(10, 4), outflow=(11, 8), where=0.45)#, patchlabel='IsoT')
     
    # IsoH in
    connection(sankey,  inflow=(10, 5), outflow=(11, 9), where=0.35)#, patchlabel='IsoH')
    
    # Dipu in
    connection(sankey, inflow=(2, 4), outflow=(11, 6), where=0.65)#, patchlabel='Dipu')
    
    sankey.finish()
    
    reservoirs = ['Räterichsboden', 'Gelmer', 'Grimsel', 'Oberaar']
    patches = ['//', '\\\\', '//', '\\\\']
    for i0, k0 in enumerate(reservoirs):
        sankey.diagrams[names.index(k0)].patch.set_hatch(patches[i0])
    
    plt.xticks(np.arange(-99, 100, gridStep*scale))
    plt.yticks(np.arange(-99, 100, gridStep*scale))
    ax.set_aspect('equal', adjustable='box', anchor='SW')
 
    if grid:
        ax.xaxis.grid(color=[0.9, 0.9, 0.9], linestyle='solid')
        ax.yaxis.grid(color=[0.9, 0.9, 0.9], linestyle='solid')
    ax.yaxis.set_ticklabels([])
    ax.xaxis.set_ticklabels([])
     
    ax.set_xlim([-10.5, 16.5])
    ax.set_ylim([-3.5, 20])
    ax.set_ylabel(title)
    
    if legend:
        leg = plt.legend(loc='center right', bbox_to_anchor=(1.55, 0.5), frameon=True, framealpha=1, fancybox=False, prop={'weight': 'normal'})
        leg.get_frame().set_linewidth(1.25)
        leg.get_frame().set_edgecolor('black')
    plt.show(block=block)
    
if __name__=='__main__':
    sns.set_palette('deep')
    sns.set(style="white", context="paper")
    
    data = pd.read_csv('KWO averages.csv', index_col=[0], header=None).iloc[:,0]
    KWOSankey(data.transpose()*365.25, title='KWO flows', palette='rainbow', details=True, legend=True, scale=1/365.25, units='Mm³/yr')
    
    plt.show(block=True)
    