#! /usr/bin/env python3
import argparse
import re,sys,os,math,gc
import numpy as np
import pandas as pd
import matplotlib as mpl
import copy
import math
from math import pi
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import sparse
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import copy
import math
import seaborn as sns
#from scipy.interpolate import BSpline, make_interp_spline

plt.rcParams.update({'figure.max_open_warning': 100000})
plt.style.use('seaborn-colorblind')
mpl.rcParams['ytick.direction'] = 'out'

mpl.rcParams['savefig.dpi'] = 300 #图片像素
mpl.rcParams['figure.dpi'] = 300 
mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']=42

__author__ ='赵玥'
__mail__   ='yuezhao@genome.cn'
_data__   ='20191101'

def draw_boundaries(ax,Boundary_dict,start,end,samplelist,str_x,sam_x):
    ax.tick_params(top='off',bottom='off',left='on',right='off')
    for loc in ['top','left','right','bottom']:
        ax.spines[loc].set_visible(False)
    #ax.spines['left'].set_color('k')
    #ax.spines['left'].set_linewidth(2)
    #ax.spines['left'].set_smart_bounds(True)
    #ax.spines['left'].set_linewidth(1)
    #ax.spines['right'].set_visible(False)
    #ax.spines['bottom'].set_visible(False)
    ax.set_axis_bgcolor('w')
    ax.set(xticks=[])
    ax.set(yticks=[])
    sample1 = samplelist[0]
    sample2 = samplelist[1]
    boundary_mid1 = Boundary_dict[sample1]['mid'].tolist()
    boundary_mid2 = Boundary_dict[sample2]['mid'].tolist()
    bound_y1min   = [1.25 for i in boundary_mid1]
    bound_y1max   = [1.75 for i in boundary_mid1]
    bound_y2min   = [0.25 for i in boundary_mid2]
    bound_y2max   = [0.75 for i in boundary_mid2]
    ax.set_ylim(0,2)
    ax.vlines(boundary_mid1,bound_y1min,bound_y1max,lw=2,color='red')
    ax.vlines(boundary_mid2,bound_y2min,bound_y2max,lw=2,color='green')
    ax.set_xlim(start,end)
    ax.text(str_x,0.5,'bound',horizontalalignment='right',verticalalignment='center',rotation='vertical',transform=ax.transAxes,fontsize=8)
    ax.text(sam_x,0.75,sample1,horizontalalignment='right',verticalalignment='center',rotation='horizontal',transform=ax.transAxes,color="red",fontsize=8)
    ax.text(sam_x,0.25,sample2,horizontalalignment='right',verticalalignment='center',rotation='horizontal',transform=ax.transAxes,color="green",fontsize=8)

def cut_boundaries(Boundary_dict,sample,boundaryPath,chrom,start,end):
    Boundary_df        = pd.read_table(boundaryPath,header=0,index_col=None,encoding='utf-8')
    Boundary_df        = Boundary_df.fillna(0)
    Boundary_df        = Boundary_df[['start','end']]
    Boundary_df['mid'] = (Boundary_df['start'] + Boundary_df['end'])/2
    Boundary_df        = Boundary_df[Boundary_df['mid']>=start]
    Boundary_df        = Boundary_df[Boundary_df['mid']<=end]
    Boundary_df.reset_index(drop=True)
    Boundary_dict[sample] = Boundary_df
    return Boundary_dict


def draw_insulation(ax,insu,chrs,start,end,color):
    #df_insu=cut_insulation(insu,chrs,start,end)
    df_insu=pd.read_table(insu,sep='\t',names=['chrs','start','end','insu'])
    ax.tick_params(top='off',bottom='off',left='on',right='off')
    line=ax.plot(df_insu['start'],df_insu['insu'], color=color, linewidth=0.8, label="insulation")
    ax.set_xlim(start,end)
    ax.set_xticks([])
    ax.set_ylim(df_insu['insu'].min(),df_insu['insu'].max())
    #ax.set_yticks([df_insu['insu'].min(),df_insu['insu'].max()])
    for loc in ['left','top','bottom']:
        ax.spines[loc].set_linewidth(0)
        ax.spines[loc].set_color('black')
    ax.spines['right'].set_linewidth(0)
    ax.spines[loc].set_color('black')

def draw_SV(files,ax,chrom,start,end,sample,color,types):
    markdf=pd.read_table(files,sep='\t')
    markdf=markdf[markdf['types']==types]
    markdf=markdf[markdf['chrs']==chrom]
    markdf=markdf[markdf['start']>start]
    markdf=markdf[markdf['end']<end]
    ax.tick_params(left='on',right='off',top='off',bottom='on')
    markdf['width'] = markdf['end'] - markdf['start']
    markdf['sign']=[1]*len(markdf)
    #vectorf = np.vectorize(np.float)
    #vectori = np.vectorize(np.int)
    #starts=list(markdf['start'])
    #hight=list(markdf['sign'])
    #width=(markdf['width'])
    ax.bar(x=list(markdf['start']),height=list(markdf['sign']),bottom=0, width = list(markdf['width']),color=color,linewidth=0,align='edge')
    ax.set_xlim([start,end])
    ax.set_ylim([0,1])
    xts   = np.linspace(start,end,2)
    yts   = np.linspace(0,1,2)
    xtkls = ['{:,}'.format(int(i)) for i in xts]
    ytkls = ['{:,}'.format(int(j)) for j in yts]
    ax.tick_params(direction='out',pad=1)
    ax.set_yticks([])
    #ax.set_yticklabels(ytkls,fontsize=5)
    ax.text(-0.11,0.0,sample,fontsize=12,color='k',horizontalalignment='left',verticalalignment='bottom',rotation='vertical',transform=ax.transAxes)
    #ax.set_title("{}_{}_{}_{}".format(sample,chrom,start,end),fontsize=10)
    ax.spines['bottom'].set_linewidth(0)
    ax.spines['left'].set_linewidth(0)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)
    if type =='bottom':
        ax.set_xticks(xts)
        ax.set_xticklabels(xtkls,fontsize=12)
        ax.spines['bottom'].set_linewidth(0.5)
        ax.spines['bottom'].set_color('k')
        ax.text(-0.11,-0.7,chrom,fontsize=12,color='k',horizontalalignment='left',verticalalignment='bottom',rotation='horizontal',transform=ax.transAxes)
    else:
        ax.set_xticks([])
        ax.set_xticklabels('')
    markdf = pd.DataFrame()
    gc.collect()

def cut_insulation(insu,chrs,start,end):
    file=open(insu)
    file_list=[]
    for i in file: 
        i=i.strip()
        file_list.append(i)
    insu_list=[]
    for i in range(len(file_list)):
        x=file_list[i].split('/')
        insu_list.append([x[-2],file_list[i]])
    list_df=pd.DataFrame(insu_list,columns=['chrs','insu'])
    list_df=list_df[list_df['chrs']==chrs]
    list_df=list_df.reset_index(drop=True)
    df_insu=pd.read_table(list_df['insu'][0],sep='\t',names=['chrs','start','end','insu'],comment='t')
    df_insu['mid']=(df_insu['start']+df_insu['end'])/2
    df_insu=df_insu.fillna(0)
    df_insu=df_insu[(df_insu['start']>start)&(df_insu['end']<end)]
    return df_insu
    


def draw_AB(files,res,chrom,start,end,sample,ax):
    compartdf = pd.read_table(files,sep='\t',names=['chrom','start','end','eigen1'])
    compartdf = compartdf[compartdf['chrom']==chrom]
    compartdf = compartdf.reset_index(drop=True)
    df = compartdf
    df=df[df['end']>=start]
    df=df[df['start']<=end]
    df=df.reset_index(drop=True)
    ax.tick_params(top='off',bottom='on',left='off',right='off')
    for loc in ['left','right','top','bottom']:
        ax.spines[loc].set_visible(False)
    df['width']=df['end']-df['start']
    #ax.axis([start, end, min,max])
    for i in range(len(df)):
        if df['eigen1'][i]>0:
            ax.bar(x=df['start'][i],height=df['eigen1'][i],bottom=0, width = df['width'][i],color='#E7605B',linewidth=0,align='edge')
        else:
            ax.bar(x=df['start'][i],height=df['eigen1'][i],bottom=0, width = df['width'][i],color='#3B679E',linewidth=0,align='edge')
    ax.set_ylim(-0.1,0.1)
    ax.set_ylabel(sample)
    ax.set_yticks([])
    ax.set_xticks([])
       

def Express_Swith(Epipath,chrom,start,end):
    Expressdf = pd.read_table(Epipath,header=None,index_col=False,sep='\t')
    Expressdf.columns = ['chrom','start','end','sign']
    Expressdf = Expressdf[Expressdf['chrom']==chrom]
    Expressdf = Expressdf[Expressdf['start']>=int(start)]
    Expressdf = Expressdf[Expressdf['end']<=int(end)]
    Expressdf = Expressdf.reset_index(drop=True)
    return Expressdf

def draw_epigenetic(file,ax,chrom,start,end,sample,color,MaxYlim,type,mins):
    markdf=pd.read_table(file,sep='\t',names=['chrs','start','end','sign'])
    markdf=markdf[markdf['chrs']==chrom]
    markdf=markdf[markdf['start']>start]
    markdf=markdf[markdf['end']<end]
    ax.tick_params(left='on',right='off',top='off',bottom='on')
    markdf['width'] = markdf['end'] - markdf['start']
    
    recs  = ax.bar(x=list(markdf['start']),height=list(markdf['sign']),bottom=0, width = list(markdf['width']),color=color,linewidth=0,align='edge')
    if MaxYlim == 'None':
        ymaxlim  = markdf['sign'].max()
        yminlim  = markdf['sign'].min()
    else:
        ymaxlim = float(MaxYlim)
        yminlim = float(mins)
    ax.set_xlim([start,end])
    ax.set_ylim([yminlim,ymaxlim])
    xts   = np.linspace(start,end,5)
    yts   = np.linspace(yminlim,ymaxlim,2)
    xtkls = ['{:,}'.format(int(i)) for i in xts]
    ytkls = ['{:,}'.format(float(j)) for j in yts]
    ax.tick_params(direction='out',pad=1)
    ax.set_yticks(yts)
    ax.set_yticklabels(ytkls,fontsize=5)
    ax.text(-0.11,0.4,sample,fontsize=6,color='k',horizontalalignment='left',verticalalignment='bottom',rotation='horizontal',transform=ax.transAxes)
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)
    #ax.set_title("{}_{}_{}_{}".format(sample,chrom,start,end),fontsize=10)
    if type =='bottom':
        ax.set_xticks(xts)
        ax.set_xticklabels(xtkls,fontsize=8)
        ax.spines['bottom'].set_linewidth(0.5)
        ax.spines['bottom'].set_color('k')
        ax.text(-0.11,-0.7,chrom,fontsize=8,color='k',horizontalalignment='left',verticalalignment='bottom',rotation='horizontal',transform=ax.transAxes)
    else:
        ax.set_xticks([])
        ax.set_xticklabels('')
    markdf = pd.DataFrame()
    gc.collect()



def draw_epigenetic2(file,ax,chrom,start,end,sample,color,MaxYlim,type,mins):
    markdf=pd.read_table(file,sep='\t',names=['chrs','start','end','sign'])
    #print (markdf.head())
    markdf=markdf[markdf['chrs']==chrom]
    markdf=markdf[markdf['start']>start]
    markdf=markdf[markdf['end']<end]
    ax.tick_params(left='on',right='off',top='off',bottom='on')
    markdf['width'] = markdf['end'] - markdf['start']
    markdf['width'] = markdf['end'] - markdf['start']
    x = np.linspace(start,end,int(len(markdf)/8))
    a_BSpline=make_interp_spline(markdf['start'],markdf['sign'],k=3)
    y_new=a_BSpline(x)
    ax.plot(x, y_new, color=color,linewidth=2)
    ax.fill_between(x,y_new ,0,facecolor=color,linewidth=0,label=sample)
    if MaxYlim == 'None':
        ymaxlim  = markdf['sign'].max()
        yminlim  = markdf['sign'].min()
    else:
        ymaxlim = float(MaxYlim)
        yminlim = float(mins)
    ax.set_xlim([start,end])
    ax.set_ylim([yminlim,ymaxlim])
    xts   = np.linspace(start,end,4)
    yts   = np.linspace(yminlim,ymaxlim,2)
    xtkls = ['{:,}'.format(int(i)) for i in xts]
    ytkls = ['{:,}'.format(int(j)) for j in yts]
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)
    ax.tick_params(top=False,right=False,width=1,colors='black',direction='out')
    ax.set_yticks(yts)
    ax.set_yticklabels(ytkls,fontsize=12)
    ax.text(-0.11,0.0,sample,fontsize=12,color='k',horizontalalignment='left',verticalalignment='bottom',rotation='vertical',transform=ax.transAxes)
    #ax.set_title("{}_{}_{}_{}".format(sample,chrom,start,end),fontsize=10)
    if type =='bottom':
        ax.set_xticks(xts)
        ax.set_xticklabels(xtkls,fontsize=12)
        ax.spines['bottom'].set_linewidth(0.5)
        ax.spines['bottom'].set_color('k')
        ax.text(-0.11,-0.7,chrom,fontsize=8,color='k',horizontalalignment='left',verticalalignment='bottom',rotation='horizontal',transform=ax.transAxes)
    else:
        ax.set_xticks([])
        ax.set_xticklabels('')
    markdf = pd.DataFrame()
    gc.collect()

def draw_RNA(file,ax,chrom,start,end,sample,color,MaxYlim,type,mins):
    markdf=pd.read_table(file,sep='\t',names=['chrs','start','end','sign'])
    #print (markdf.head())
    markdf=markdf[markdf['chrs']==chrom]
    markdf=markdf[markdf['start']>start]
    markdf=markdf[markdf['end']<end]
    ax.tick_params(left='on',right='off',top='off',bottom='on')
    markdf['width'] = markdf['end'] - markdf['start']
    vectorf = np.vectorize(np.float)
    vectori = np.vectorize(np.int)
    starts=vectori(markdf['start'])
    hight=vectorf(markdf['sign'])
    width=vectori(markdf['width'])
    ax.bar(x=starts,height=hight,bottom=0,width=width,color=color,linewidth=0,align='edge')
    if MaxYlim == 'None':
        ymaxlim  = markdf['sign'].max()
        yminlim  = markdf['sign'].min()
    else:
        ymaxlim = float(MaxYlim)
        yminlim = float(mins)
    ax.set_xlim([start,end])
    ax.set_ylim([yminlim,ymaxlim])
    xts   = np.linspace(start,end,5)
    yts   = np.linspace(yminlim,ymaxlim,2)
    xtkls = ['{:,}'.format(int(i)) for i in xts]
    ytkls = ['{:,}'.format(int(j)) for j in yts]
    ax.tick_params(direction='out',pad=1)
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)
    ax.set_yticks(yts)
    ax.set_yticklabels(ytkls,fontsize=12)
    ax.text(-0.11,0.4,sample,fontsize=12,color='k',horizontalalignment='left',verticalalignment='bottom',rotation='vertical',transform=ax.transAxes)
    #ax.set_title("{}_{}_{}_{}".format(sample,chrom,start,end),fontsize=10)
    if type =='bottom':
        ax.set_xticks(xts)
        ax.set_xticklabels(xtkls,fontsize=12)
        ax.spines['bottom'].set_linewidth(0.5)
        ax.spines['bottom'].set_color('k')
        ax.text(-0.11,-0.7,chrom,fontsize=12,color='k',horizontalalignment='left',verticalalignment='bottom',rotation='horizontal',transform=ax.transAxes)
    else:
        ax.set_xticks([])
        ax.set_xticklabels('')
    markdf = pd.DataFrame()
    gc.collect()



def Express_Swith(Epipath,chrs,start,end):
    Expressdf = pd.read_table(Epipath,header=None,index_col=False,sep='\t')
    Expressdf.columns = ['chrs','start','end','sign']
    Expressdf = Expressdf[Expressdf['chrs']==chrs]
    Expressdf = Expressdf[Expressdf['start']>=int(start)]
    Expressdf = Expressdf[Expressdf['end']<=int(end)]
    Expressdf = Expressdf.reset_index(drop=True)
    return Expressdf

def draw_diff_epigenetic(file1,file2,ax,chrs,start,end,color,MaxYlim,MinYlim,type):
    df1=Express_Swith(file1,chrs,start,end)
    df2=Express_Swith(file2,chrs,start,end)
    markdf = pd.merge(df1,df2,on='start',how='inner')
    markdf['sign'] = np.log2(markdf['sign_x']) - np.log2(markdf['sign_y'])
    markdf         = markdf[['chrs_x','start','end_x','sign']]
    markdf.columns = ['chrs','start','end','sign']
    markdf         = markdf.reset_index(drop=True)
    ax.tick_params(left='on',right='off',top='off',bottom='on')
    markdf['width'] = markdf['end'] - markdf['start']
    recs  = ax.bar(markdf['start'],markdf['sign'],bottom=0, width = markdf['width'],color=color,linewidth=0)
    if MaxYlim == 'None':
        ymaxlim  = markdf['sign'].max()
        yminlim  = markdf['sign'].min()
    else:
        ymaxlim = float(MaxYlim)
        yminlim = float(MinYlim)
    ax.set_xlim([start,end])
    ax.set_ylim([yminlim,ymaxlim])
    xts   = np.linspace(start,end,5)
    yts   = np.linspace(yminlim,ymaxlim,2)
    xtkls = ['{:,}'.format(int(i)) for i in xts]
    ytkls = ['{:,}'.format(int(j)) for j in yts]
    
    ax.tick_params(direction='out',pad=1)
    ax.set_yticks(yts)
    ax.set_yticklabels(ytkls,fontsize=5)
    #ax.text(-0.11,0.4,sample,fontsize=6,color='k',horizontalalignment='left',verticalalignment='bottom',rotation='horizontal',transform=ax.transAxes)
    #ax.set_title("{}_{}_{}_{}".format(sample,chrom,start,end),fontsize=10)
    if type =='bottom':
        ax.set_xticks(xts)
        ax.set_xticklabels(xtkls,fontsize=8)
        ax.spines['bottom'].set_linewidth(0.5)
        ax.spines['bottom'].set_color('k')
        ax.text(-0.11,-0.7,chrs,fontsize=8,color='k',horizontalalignment='left',verticalalignment='bottom',rotation='horizontal',transform=ax.transAxes)
    else:
        ax.set_xticks([])
        ax.set_xticklabels('')
    markdf = pd.DataFrame()
    ax.spines['bottom'].set_linewidth(0)
    ax.spines['left'].set_linewidth(1)
    ax.spines['top'].set_linewidth(0)
    ax.spines['right'].set_linewidth(0)
    gc.collect()






def draw_bar(ax,file,chrom,start,end,max,min):
    df=pd.read_table(file,sep='\t',names=['chrs','start','end','sign'])
    df=df[df['chrs']==chrom]
    df=df[df['start']>start]
    df=df[df['end']<end]
    df=df.reset_index(drop=True)
    ax.tick_params(top='off',bottom='on',left='off',right='off')
    for loc in ['left','right','top']:
        ax.spines[loc].set_visible(False)
    df['width']=df['end']-df['start']
    #ax.axis([start, end, min,max])
    for i in range(len(df)):
        if df['sign'][i]>0:
            ax.bar(df['start'][i],df['sign'][i],bottom=0, width = df['width'][i],color='#E7605B',linewidth=0)
        else:
            ax.bar(df['start'][i],df['sign'][i],bottom=0, width = df['width'][i],color='#3B679E',linewidth=0)
    ax.set_ylim(min,max)
    ax.set_yticks([])
    ax.set_xticks([])

def get_4C_data(matrix,tstart,tend,binsize,start,end):
    print (binsize)
    t=int((tstart-start)/int(binsize))
    print ('t',t,'matrix',len(matrix))
    datalist=matrix.loc[:,[t]]
    return datalist

from statsmodels.nonparametric.smoothers_lowess import lowess
def draw_4C_module(ax,df_list,chrs,start,end,color_list,ymin,ymax,sample_list):
    ax.tick_params(top='off',bottom='off',left='on',right='off')
    i=0 
    for df in df_list:
        x = np.linspace(start,end,len(df))
        df['width']=df['end']-df['start']
        df_loess = pd.DataFrame(lowess(df['sign'], np.arange(len(df['sign'])), frac=0.05)[:, 1], index=df.index, columns=['sign'])
        ax.plot(x,df_loess['sign'], color=color_list[i], linewidth=2,label=sample_list[i],alpha=0.3)
        i+=1
        #ax.legend(handles2, labels2)
    ax.set_xlim(start,end)
    ax.set_ylim(ymin,ymax)
    ax.set_yticks([ymin,ymax])
    ax.legend(loc='right',bbox_to_anchor=(1.05,0.3),handlelength=1,handleheight=0.618,fontsize=6,frameon=False)
    for loc in ['left']:
        ax.spines[loc].set_linewidth(0.6)
        ax.spines[loc].set_color('gray')
    #ax.tick_params(top=False,right=False,width=1,colors='black',direction='out')
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)
    ax.tick_params(top=False,right=False,bottom=False,width=1,colors='black',direction='out')

def draw_4C(ax,chrs,start,end,matrix_list,samples,binsize,tstart,tend,colors,ymax):
    sample_list=samples.split(',')
    bed_list=[]
    for sample in sample_list:
        matrix,min=extract_raw_matrix(matrix_list,sample,chrs,start,end,binsize)
        datalist=get_4C_data(matrix,int(tstart),int(tend),binsize,int(start),int(end))
        bed_list.append(datalist)
    starts=[]
    for i in range(start,end,int(binsize)):
        starts.append(i)
    df_list=[]
    for i in bed_list:
        df=pd.DataFrame({'start':starts})
        df['chrs']=[chrs]*len(df)
        df['end']=df['start']+int(binsize)
        df['sign']=i
        df_list.append(df)
    color_list=colors.split(',')
    draw_4C_module(ax,df_list,chrs,start,end,color_list,0,int(ymax),sample_list)

def draw_compartment(ax,sample,compmergedf,chrom,start,end,type='top'):
    ax.tick_params(top='off',bottom='on',left='off',right='off')
    for loc in ['left','right','top']:
        ax.spines[loc].set_visible(False)
    mat    = compmergedf[sample]
    #print(mat)
    s      = compmergedf['start']
    colors = ['red','blue','#458B00','#B9BBF9','black']
    ax.set_xlim(start,end)
    if sample == 'Merge':
        ax.fill_between(s, 0, 0.25,where=mat==1, facecolor=colors[0],linewidth=0,label='CompartmentA')
        ax.fill_between(s, 0.25, 0.5,where=mat==2, facecolor=colors[2],linewidth=0,label='A Switch B')
        ax.fill_between(s, 0, -0.25,where=mat==-1,facecolor=colors[1],linewidth=0,label='CompartmentB')
        ax.fill_between(s, -0.25,-0.5,where=mat==-2,facecolor=colors[3],linewidth=0,label='B Switch A')
        legend = ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.,prop={'size':4},ncol=1)
        legend.get_frame().set_facecolor('white')
    else:
        ax.fill_between(s, 0, mat,where=mat>= 0, facecolor=colors[0],linewidth=0,label='CompartmentA')
        ax.fill_between(s, 0, mat,where=mat< 0, facecolor=colors[1],linewidth=0,label='CompartmentB')
    #ax.text(max(mat)/4,-5,'A');ax.text(max(mat)/2,-5,'B')
    ax.text(-0.11,0.4,sample,fontsize=6,color='k',horizontalalignment='left',verticalalignment='bottom',rotation='horizontal',transform=ax.transAxes)
    ymax  = mat.max()+0.005
    ymin  = mat.min()-0.005
    xts   = np.linspace(start,end,5)
    xtkls = ['{:,}'.format(int(i)) for i in xts]
    ax.set_ylim(ymin,ymax)
    ax.set_yticks([])
    #ax.set_ylabel(sample,rotation='vertical',fontsize='small')
    #compmergedf = pd.DataFrame()
    if type =='bottom':
        ax.set_xticks(xts)
        ax.set_xticklabels(xtkls,fontsize=8)
        ax.spines['bottom'].set_linewidth(1)
        ax.spines['bottom'].set_color('k')
        ax.text(-0.11,-0.7,chrom,fontsize=8,color='k',horizontalalignment='left',verticalalignment='bottom',rotation='horizontal',transform=ax.transAxes)
    else:
        ax.set_xticks([])
        ax.spines['bottom'].set_visible(False)
    gc.collect()



def colorbar(ax,im,vmin,vmax):
    axins1 = inset_axes(ax, width=0.1,height=0.6,loc=3, bbox_to_anchor=(0, 0.2, 0.5, 1), bbox_transform=ax.transAxes,borderpad=0)
    print (vmin,vmax)
    cbar=plt.colorbar(im, cax=axins1, orientation='vertical',ticks=[math.ceil(vmin),int(vmax)])
    axins1.tick_params(left='on',right='off',top='off',bottom='off',labelsize=12)
    axins1.yaxis.set_ticks_position('left')
    return cbar

import math
from matplotlib import pyplot as plt
plt.style.use('seaborn-colorblind')
pd.set_option('display.precision',2)
from scipy import sparse
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#from annotation_GenePred import Gene_annotation


def cut_mat(mat,start,end,resolution,min):
    start = int(int(start)/resolution)
    end = math.ceil(int(end)/resolution)
    start = int(start - min)
    end = int(end - min)
    mat = mat.fillna(0)
    mat = mat.iloc[start:end+1,start:end+1]
    gc.collect()

    return mat,start,end

def self_zscore(df):
    dsc = pd.DataFrame(np.ravel(df)).describe(include=[np.number])
    df = (df - dsc.loc['mean',0])/dsc.loc['std',0]
    return df

from scipy.ndimage import gaussian_filter
def get_matrix(mat_path,binsize,start,end):
    binsize=int(binsize)
    mat=pd.read_table(mat_path,names=['b1','b2','contacts'])
    mat=mat[(mat['b1']>=start-3000000) & (mat['b2']>=start-3000000)]
    mat=mat[(mat['b1']<=end+3000000) & (mat['b2']<=end+3000000)]
    #-----------xlim genome start genome end-------------------------------
    min=mat['b1'].min()
    max=mat['b1'].max()
    min=math.ceil(int(min)/binsize)*binsize
    max=int(int(max)/binsize)*binsize
    N=int(max/binsize)-math.ceil(min/binsize)+1
    mat['b1']=mat['b1'].apply(lambda x: (x-min-1)/binsize)
    mat['b2']=mat['b2'].apply(lambda x: (x-min-1)/binsize)
    #-----------coo matrix-----------------------------------------------
    counts=sparse.coo_matrix((mat['contacts'],(mat['b1'],mat['b2'])),shape=(N, N),dtype=float).toarray()
    diag_matrix=np.diag(np.diag(counts))
    counts=counts.T + counts
    #counts=counts-diag_matrix
    counts=counts-diag_matrix-diag_matrix
    df=pd.DataFrame(counts)
    #----------zscore minus ---------------------------------
    df=self_zscore(df)
    min=int(min/binsize)
    df,min,max=cut_mat(df,start,end,binsize,min)
    np.fill_diagonal(df.values, 0)
    return df,min
    
def get_matrix_df(lists,sample,chrs):
    df=pd.read_table(lists,sep='\t',names=['sample','chrs','matrix'])
    df=df[df['sample']==sample]
    df=df[df['chrs']==chrs]
    df=df.reset_index(drop=True)
    matrix=df['matrix'][0]
    return matrix

def extract_matrix(lists,sample1,sample2,chrs,start,end,binsize):
    matrix1=get_matrix_df(lists,sample1,chrs)
    matrix2=get_matrix_df(lists,sample2,chrs)
    zscore1,min=get_matrix(matrix1,binsize,start,end)
    zscore2,min=get_matrix(matrix2,binsize,start,end)
    delta=zscore1-zscore2
    delta = gaussian_filter(delta, sigma=0.5)
    #start=(start - min)/binsize
    #end= (end -min)/binsize
    #delta=delta.loc[start:end,start:end]
    return delta,min


class MidpointNormalize(mpl.colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)
    
    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

def heat_map(ax,matrix,title,start,end):
    for loc in ['left','right','top']:
        ax.spines[loc].set_visible(False)
    matrix = Triangular(matrix)
    triArr = np.array(matrix)
    triArr_ravel = pd.DataFrame(triArr.ravel())
    triArr_ravel = triArr_ravel[triArr_ravel!= -np.inf]
    triArr_ravel = triArr_ravel[triArr_ravel!=np.inf]
    #vmax=np.nanpercentile(triArr_ravel,97)
    #vmin=np.nanpercentile(triArr_ravel,0)
    triMat       = ''
    triArr_ravel = ''
    tmp_trimat   = ''
    matrix= triArr
    my_colors=['#5A1216','#A61B29','#F0A1A8','#E3B4B8','#FFFEF8','#93B5CF','#2775B6','#144A74','#101F30']
    colormap = mpl.colors.LinearSegmentedColormap.from_list('cmap',my_colors[::-1],500)
    pp = matrix[matrix>0]
    pm = matrix[matrix<0]
    vmax = np.nanpercentile(pp, 99)
    vmin = np.nanpercentile(pp, 1)
    im=ax.imshow(matrix,cmap=colormap,interpolation="nearest",aspect='auto',norm=MidpointNormalize(midpoint=0,vmin=-1, vmax=1),origin='lower')
    colormap.set_bad('white')
    ax.set_title(title,y=0.5,x=0,fontsize=12)
    ax.tick_params(direction='out',pad=5)
    ax.spines['bottom'].set_linewidth(0)
    ax.spines['left'].set_linewidth(0)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)
    ax.tick_params(top=False,right=False,left=False,bottom=False,width=0,colors='black',direction='out')
    ax.set_yticks([])
    n=len(matrix)
    #ax.set_xticks([0,n*2])
    ax.set_xticks([])
    ax.set_xticklabels([str(start),str(end)])
    cbar=colorbar(ax,im,vmin,vmax)


def get_raw_matrix(mat_path,binsize,start,end):
    mat=pd.read_table(mat_path,names=['frag1','frag2','contacts'])
    mat=mat[(mat['frag1']>=start) & (mat['frag2']>=start)]
    mat=mat[(mat['frag1']<=end) & (mat['frag2']<=end)]
    #-----------xlim genome start genome end-------------------------------
    min=mat['frag1'].min()
    max=mat['frag1'].max()
    binsize=int(binsize)
    min=math.ceil(int(min)/binsize)*binsize
    max=int(int(max)/binsize)*binsize
    N=int(max/binsize)-math.ceil(min/binsize)+1
    print (N)    

    #------------------tranform matrix ----------------------------------
    mat['b1']=mat['frag1'].apply(lambda x: (x-min)/binsize)
    mat['b2']=mat['frag2'].apply(lambda x: (x-min)/binsize)
    #-----------coo matrix-----------------------------------------------
    counts=sparse.coo_matrix((mat['contacts'],(mat['b1'],mat['b2'])),shape=(N, N),dtype=float).toarray()
    diag_matrix=np.diag(np.diag(counts))
    counts=counts.T + counts
    #counts=counts-diag_matrix
    counts=counts-diag_matrix
    df=pd.DataFrame(counts)
    #----------zscore minus ---------------------------------
    np.fill_diagonal(df.values, 0)
    return df,min

def extract_raw_matrix(lists,sample,chrs,start,end,binsize):
    matrix=get_matrix_df(lists,sample,chrs)
    matrix,min=get_raw_matrix(matrix,binsize,start,end)
    return matrix,min

def Triangular(mat):
    mat=np.array(mat)
    n=0
    length=len(mat)
    tri_mat=np.zeros([length,length*2])
    tri_mat[tri_mat==0]=np.nan
    for i in range(length):
        curl=np.array(np.diag(mat,i))
        tri_mat[i,n:(length*2-n)]=curl.repeat(2)
        n+=1
    mat  = ''
    gc.collect()
    return tri_mat

def draw_delta(ax,chrs,start,end,matrix_list,samples,binsize):
    samples=samples.split(',')
    sample1,sample2=samples[0],samples[1]
    delta,min=extract_matrix(matrix_list,sample1,sample2,chrs,start,end,binsize)
    heat_map(ax,delta,'{}-\n{}'.format(sample1,sample2),start,end)
    #return datalist

def rawheat_map(ax,matrix,title,start,end):
    #fig = plt.figure(figsize=(6,6))
    #ax = fig.add_axes([0.15,0.15,0.7,0.7])
    for loc in ['left','right','top']:
        ax.spines[loc].set_visible(False)
    matrix = Triangular(matrix)
    triArr = np.array(matrix)
    triArr_ravel = pd.DataFrame(triArr.ravel())
    triArr_ravel = triArr_ravel[triArr_ravel!= -np.inf]
    triArr_ravel = triArr_ravel[triArr_ravel!=np.inf]
    vmax=np.nanpercentile(triArr_ravel,97)
    vmin=np.nanpercentile(triArr_ravel,0)
    triMat       = ''
    triArr_ravel = ''
    tmp_trimat   = ''
    matrix= triArr
    my_colors=['#5A1216','#C02C38','#EE3F4D','#F07C82','#F1FB00','#2AFFBD','#0000A4']
    colormap = mpl.colors.LinearSegmentedColormap.from_list('cmap',my_colors[::-1],500)
    colormap.set_bad('white')
    im=ax.imshow(matrix,cmap=colormap,clim=(vmin,vmax),interpolation="nearest",aspect='auto',origin='lower')
    ax.set_title(title,y=0.5,x=0,fontsize=16)
    ax.tick_params(direction='out',pad=5)
    ax.spines['bottom'].set_linewidth(0)
    ax.spines['left'].set_linewidth(0)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)
    ax.tick_params(top=False,right=False,left=False,bottom=False,width=0,colors='black',direction='out')
    #ticks=np.linspace(0,length*2,3)
    n=len(matrix)
    #ax.set_xticks([0,n*2])
    ax.set_xticks([])
    #ax.set_xticklabels([str(start),str(end)])
    ax.set_yticks([])
    #ax.set_xticklabels(sampleLst)
    #vmin=matrix.stack().min()
    #vmax=matrix.stack().max()
    print (ax)
    cbar=colorbar(ax,im,vmin,vmax)
    #fig.savefig('{}.pdf'.format(outfile))

def draw_matrix(ax,chrs,start,end,matrix_list,sample,binsize):
    matrix,min=extract_raw_matrix(matrix_list,sample,chrs,start,end,binsize)
    matrix=gaussian_filter(matrix, sigma=0.5)
    rawheat_map(ax,matrix,sample,start,end)
    #return

def draw_compare_border(sample,ax,files,chrs,start,end,color):
    df=pd.read_table(files,sep='\t',names=['chrs','start','end'])
    df=df[df['chrs']==chrs]
    df=df[df['start']>start]
    df=df[df['end']<end]
    df['sign']=1
    df['width']=(end-start)/198
    title = sample
    ax.axis([start, end, 0,1])
    ax.bar(x=list(df['start']),height=list(df['sign']),bottom=0,width=list(df['width']),color=color,align='edge')
    #ax.bar(df['start'],df['sign'],bottom=0, width = df['width'],color=color,linewidth=0)
    for loc in ['bottom','left','right','top']:
        ax.spines[loc].set_linewidth(0)  
    ax.set_xticks([])
    ax.set_yticks([])

def draw_CDB(sample,ax,files,chrs,start,end):
    df=pd.read_table(files,sep='\t')
    df=df[df['chrs']==chrs]
    df=df[df['start']>start]
    df=df[df['end']<end]
    df['sign']=[1]*len(df)
    df['width']=(end-start)/198
    df1=df[df['level']==0]
    df2=df[df['level']==1]
    title = sample
    df1=df1.reset_index(drop=True)
    df2=df2.reset_index(drop=True)
    df1=df1.loc[:,['chrs','start','end','sign','width','level']]
    df2=df2.loc[:,['chrs','start','end','sign','width','level']]
    ax.axis([start, end, 0,1])
    ax.bar(df1['start'],df1['sign'],bottom=0, width = df1['width'],color='#5A685B',linewidth=0)
    ax.bar(df2['start'],df2['sign'],bottom=0, width = df2['width'],color='#FEC79E',linewidth=0)
    ax.set_xticks([])
    ax.set_yticks([])



def gene_match(geneinfo):

    gene_version ="";gene_name="";
    #gene_id "ENSRNOG00000046319"; gene_version "4";
    #gene_name "Vom2r3"; gene_source "ensembl_havana";
    match   = re.search("gene_name \"(.+?)\";",geneinfo)
    if match:
        gene_name   = match.group(1)
        #print(gene_name)
        return (gene_name)

def transcript_match(geneinfo):
    #transcript_id "ENSRNOT00000044187";
    transcript_id = ''
    match = re.search("transcript_id \"(.+?)\";",geneinfo)
    if match:
        transcript_id   = match.group(1)
        return (transcript_id)

def exon_match(geneinfo):
    exon_num = 1
    match    = re.search("exon_number \"(.+?)\";",geneinfo)
    if match:
        exon_num = int(match.group(1))
    return (exon_num)
def gene_biotype_match(geneinfo):
    #gene_biotype "processed_transcript";
    gene_biotype =''
    match   =  re.search("gene_biotype \"(.+?)\";",geneinfo)
    if match:
        gene_biotype = match.group(1)
    return gene_biotype

def filter_gtffile(gtffile,chr,start,end,Type='all'):
    gtf_bed  = pd.read_table(gtffile,usecols=[0,2,3,4,6,8],names=['chr','transcript_ID','start','end','orient','desc'])

    gtf_bed  = gtf_bed[(gtf_bed['chr']==chr) & ((gtf_bed['transcript_ID']=='exon')|(gtf_bed['transcript_ID']=='three_prime_utr') | (gtf_bed['transcript_ID']=='five_prime_utr'))]
    gtf_bed = gtf_bed.reset_index(drop=True)
    start =int(start)
    end = int(end)
    gtf_bed['start'] = pd.to_numeric(gtf_bed.start, errors='coerce' )
    gtf_bed['end'] = pd.to_numeric(gtf_bed.end, errors='coerce' )
    gtf_bed = gtf_bed[(gtf_bed['start']>start) & (gtf_bed['end']<end)]
    gtf_bed = gtf_bed.reset_index(drop=True)
    gtf_bed['mid']             = (gtf_bed['start'] + gtf_bed['end'])/2
    gtf_bed['gene_desc']       = gtf_bed['desc'].apply(gene_match)
    gtf_bed['transcript_desc'] = gtf_bed['desc'].apply(transcript_match)
    gtf_bed['exon_num']        = gtf_bed['desc'].apply(exon_match)
    gtf_bed['gene_type']       = gtf_bed['desc'].apply(gene_biotype_match)
    if Type == 'all':
        gtf_bed = gtf_bed
    else:
        geneType= Type.split(' ')
        gtf_bed  = gtf_bed[gtf_bed['gene_type'].isin(geneType)]
        #gtf_bed.loc[(gtf_bed['transcript_desc'=='three_prime_utr'),'transcript_ID']='3UTR'
        #gtf_bed.loc[(gtf_bed['transcript_desc'=='five_prime_utr'),'transcript_ID']='3UTR'
    return gtf_bed

def draw_genes(ax,gtf_bed,start,end):
    ax.tick_params(left='off',top='off',bottom='on',right='off')
    ax.spines['bottom'].set_color('k')
    ax.spines['bottom'].set_linewidth(0.05)
    for loc in ['top','left','right']:
        ax.spines[loc].set_visible(False)
    gtf_bed  = gtf_bed.sort_values(['mid'],ascending=[True])
    #print(gtf_bed['gene_desc'].tolist())
    exon_bed = gtf_bed[gtf_bed['transcript_ID']=='exon']
    UTR3_bed = gtf_bed[gtf_bed['transcript_ID']=='three_prime_utr']
    UTR5_bed = gtf_bed[gtf_bed['transcript_ID']=='five_prime_utr']

    k        = 0
    genelist = []
    for i in gtf_bed.index:

        if (gtf_bed.loc[i,'gene_desc']!=None) and (gtf_bed.loc[i,'gene_desc'] not in genelist):
            genelist.append(gtf_bed.loc[i,'gene_desc'])
    cutbed   = pd.DataFrame()
    print (gtf_bed['desc'][0])
    for name in genelist:
        #print(name)
        tmpdf           = copy.deepcopy(gtf_bed)
        cutbed          = tmpdf[tmpdf['gene_desc']==name]

        exons           = np.array(cutbed[['start','end']])
        Sample_UTR3_bed = UTR3_bed[UTR3_bed['gene_desc']==name]
        Sample_UTR5_bed = UTR5_bed[UTR5_bed['gene_desc']==name]
        Sample_CDS_bed  = exon_bed[exon_bed['gene_desc']==name]
        utr3            = np.array(Sample_UTR3_bed[['start','end']])
        utr5            = np.array(Sample_UTR5_bed[['start','end']])
        cds             = np.array(Sample_CDS_bed[['start','end']])
        gene_start      = np.min(exons.ravel())
        gene_end        = np.max(exons.ravel())#143084

        genenum = len(genelist)
        num     = math.ceil(genenum/10)+2
        m       = math.fmod(k,num)
        k       = k + 1
        rand    = math.ceil(k/num)%2
        range = 1/num
        rec=mpatches.Rectangle((gene_start,0.9475 - m*range),abs(gene_end-gene_start),0.005,color='grey',alpha=0.5)
        ax.add_patch(rec)
        if rand == 1:
            ax.text(((gene_start+gene_end)/2),(0.99 - m*range),"{}".format(name),fontsize=9,horizontalalignment='center',verticalalignment='center',rotation='horizontal',color="black")
        else:
            ax.text(((gene_start+gene_end)/2),(0.905 - m*range),"{}".format(name),fontsize=9,horizontalalignment='center',verticalalignment='center',rotation='horizontal',color="black")
        for exon in cds:
            rec0  = mpatches.Rectangle((exon[0], 0.9375 - m*range),(exon[1]-exon[0]),0.025,color="black")
            ax.add_patch(rec0)
        for i in utr3:
            rec1  = mpatches.Rectangle((i[0],0.94375 - m*range ),(i[1]-i[0]),0.0125 ,color="grey",alpha=0.8)
            ax.add_patch(rec1)
        for j in utr5:
            rec2  = mpatches.Rectangle((j[0],0.94375 - m*range ),(j[1]-j[0]),0.0125 ,color="blue",alpha=0.8)
            ax.add_patch(rec2)

    xmin = start
    xmax = end
    ax.set_xlim(xmin,xmax)
    xtick  = np.linspace(start,end,4)
    xtlabs=[]
    for i in xtick:
        i=i/1000000
        i='%.2f' % i
        xtlabs.append(str(i))
    #xtlabs = ["{:,}".format(int(x)/1000000) for x in xtick]
    ax.set_ylim([0,1])
    ax.set_xticks(xtick)
    ax.set_xticklabels(xtlabs,fontsize=15)
    ax.set_yticks([])

def pairwise(file,ax,chrs,start,end,label,color):
    t = np.linspace(0, 1*pi, 100)
    df=pd.read_table(file,sep='\t',names=['chr','loci1','loci2','observed','expected'])
    df_new=df[df['chr']==chrs]
    df_new=df_new[df_new['loci1']>int(start)]
    df_new=df_new[df_new['loci2']<int(end)]
    max=(df_new['loci2']-df_new['loci1']).max()+100000
    ax.axis([start, end, 0,max])
    for i in range(len(df['chr'])):
        if str(df['chr'][i])==chrs and df['loci1'][i]>int(start) and df['loci2'][i]<int(end):
            u=(df['loci1'][i]+df['loci2'][i])/2
            v=0
            a=(df['loci2'][i]-df['loci1'][i])/2
            b=df['loci2'][i]-df['loci1'][i]
            ax.plot( u+a*np.cos(t) , v+b*np.sin(t) ,color=color,linewidth=1)
    ax.set_ylabel(label,fontsize=20)
    if math.isnan(max):
        max=100
    top=max
    xticks=[]
    print (start,end)
    for i in range(start,end,int((end-start)/5)):
        xticks.append(i)
    print (xticks)
    xticks.pop(0)
    yticks=[0,top]
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks,fontsize=3)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks,fontsize=3)
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)

def draw_line(file,ax,chrom,start,end,sample,color,MaxYlim,type,mins):
    markdf=pd.read_table(file,sep='\t',names=['chrs','start','end','sign'])
    markdf=markdf[markdf['chrs']==chrom]
    markdf=markdf[markdf['start']>start]
    markdf=markdf[markdf['end']<end]
    ax.tick_params(left='on',right='off',top='off',bottom='on')
    markdf['width'] = markdf['end'] - markdf['start']
    vectorf = np.vectorize(np.float)
    vectori = np.vectorize(np.int)
    starts=vectori(markdf['start'])
    hight=vectorf(markdf['sign'])
    width=vectori(markdf['width'])
    ax.bar(x=starts,height=hight,bottom=0,width=width,color=color,linewidth=0.1,align='edge',fill=False,edgecolor=color)
    if MaxYlim == 'None':
        ymaxlim  = markdf['sign'].max()
        yminlim  = markdf['sign'].min()
    else:    
        ymaxlim = float(MaxYlim)
        yminlim = float(mins)
    ax.set_xlim([start,end])
    ax.set_ylim([yminlim,ymaxlim])
    xts   = np.linspace(start,end,5)
    yts   = np.linspace(yminlim,ymaxlim,2)
    xtkls = ['{:,}'.format(int(i)) for i in xts] 
    ytkls = ['{:,}'.format(float(j)) for j in yts] 
    ax.tick_params(direction='out',pad=1)
    ax.tick_params(direction='out',pad=1)
    ax.set_yticks(yts)
    ax.set_yticklabels(ytkls,fontsize=5)
    ax.text(-0.11,0.4,sample,fontsize=6,color='k',horizontalalignment='left',verticalalignment='bottom',rotation='horizontal',transform=ax.transAxes)
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)
    ax.tick_params(top=False,right=False,width=1,colors='black',direction='out')    #ax.set_title("{}_{}_{}_{}".format(sample,chrom,start,end),fontsize=10)
    if type =='bottom':
        ax.set_xticks(xts)
        ax.set_xticklabels(xtkls,fontsize=8)
        ax.spines['bottom'].set_linewidth(0.5)
        ax.spines['bottom'].set_color('k')
        ax.text(-0.11,-0.7,chrom,fontsize=8,color='k',horizontalalignment='left',verticalalignment='bottom',rotation='horizontal',transform=ax.transAxes)
    else:
        ax.set_xticks([])
        ax.set_xticklabels('')
    markdf = pd.DataFrame()
    gc.collect()


def main():
    parser=argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
    parser.add_argument('-fl','--file_list',dest='file_list',type=str,required=True)
    parser.add_argument('-u','--uppercent',help='up_percent : 85',dest='uppercent',type=str,default='85')
    parser.add_argument('-d','--downpercent',help='down_percent :5',dest='downpercent',type=str,default='5')
    parser.add_argument('-o','--output',dest='output',type=str,required=True)
    parser.add_argument('-chr','--chrom',dest='chromosome',type=str,required=True)
    parser.add_argument('-st','--start',dest='start',type=str,required=True)
    parser.add_argument('-ed','--end',dest='end',type=str,required=True)
    parser.add_argument('-t','--genetype',help='gene type',dest='genetype',type=str,default='all')
    parser.add_argument('-b1','--bound1',help='bound1',dest='bound1',type=int,default='0')
    parser.add_argument('-b2','--bound2',help='bound2',dest='bound2',type=int,default='0')
    args = parser.parse_args()
    uppercent = int(args.uppercent)
    downpercent = int(args.downpercent)
    start=int(args.start)-args.bound1
    end=int(args.end)+args.bound2
    tstart=int(args.start)
    tend=int(args.end)
    chrs=args.chromosome
    files=open(args.file_list)
    file_list=[]
    for i in files:
        i=i.strip().split('\t')
        file_list.append(i)
    dic={'matrix':3,'delta':3,'gtf':3,'bedgraph':0.8,'pairwise':1.2,'SV':0.5,'insu':0.8,'delta_epi':1.6,'DI':0.8,'AB':0.8,'boundary':0.8,'assembly':1.3,'CDB':0.8,'RNA':0.8,'4C':1.2,'line':0.8}
    l=0
    for i in file_list:
        l+=dic[i[0]]+dic[i[0]]/3
    fig = plt.figure(figsize=(8,l))
    x=0
    j=0
    while j < len(file_list):
        #print (j)
        if file_list[j][0]=='matrix':
            x+=dic[file_list[j][0]]
            ax = fig.add_axes([0.1,1-(x/l)-x/(4*l),0.8,(1/l)*3.0])
            draw_matrix(ax,chrs,start,end,file_list[j][3],file_list[j][1],file_list[j][2])
            j+=1
        elif file_list[j][0]=='delta':
            x+=dic[file_list[j][0]]
            ax = fig.add_axes([0.1,1-(x/l)-x/(4*l),0.8,(1/l)*3.0])
            draw_delta(ax,chrs,start,end,file_list[j][3],file_list[j][1],file_list[j][2])
            j+=1
        elif file_list[j][0]=='4C':
            x+=dic[file_list[j][0]]
            ax = fig.add_axes([0.1,1-(x/l)-x/(4*l),0.8,(1/l)*1.1])
            #matrix_list:4,samples:2,colors:3,binsize:1,ymax:5
            draw_4C(ax,chrs,start,end,file_list[j][4],file_list[j][2],file_list[j][1],tstart,tend,file_list[j][3],file_list[j][5])
            j+=1
        elif file_list[j][0]=='bedgraph':
            x+=dic[file_list[j][0]]
            ax = fig.add_axes([0.1,1-(x/l)-x/(4*l),0.8,(1/l)*0.65])
            draw_epigenetic(file_list[j][3],ax,chrs,start,end,file_list[j][1],file_list[j][2],file_list[j][4],'up',file_list[j][5])
            j+=1
        elif file_list[j][0]=='line':
            x+=dic[file_list[j][0]]
            ax = fig.add_axes([0.1,1-(x/l)-x/(4*l),0.8,(1/l)*0.65])
            draw_line(file_list[j][3],ax,chrs,start,end,file_list[j][1],file_list[j][2],file_list[j][4],'up',file_list[j][5])
            j+=1 
        elif file_list[j][0]=='RNA':
            x+=dic[file_list[j][0]]
            ax = fig.add_axes([0.1,1-(x/l)-x/(4*l),0.8,(1/l)*0.65])
            draw_RNA(file_list[j][3],ax,chrs,start,end,file_list[j][1],file_list[j][2],file_list[j][4],'up',file_list[j][5])
            j+=1
        elif file_list[j][0]=='SV':
            x+=dic[file_list[j][0]]
            ax = fig.add_axes([0.1,1-(x/l)-x/(4*l),0.8,(1/l)*0.35])
            draw_SV(file_list[j][4],ax,chrs,start,end,file_list[j][1],file_list[j][2],file_list[j][3])
            j+=1
        elif file_list[j][0]=='pairwise':
            x+=dic[file_list[j][0]]
            ax = fig.add_axes([0.1,1-(x/l)-x/(4*l),0.8,(1/l)*1.1])
            pairwise(file_list[j][3],ax,chrs,start,end,file_list[j][1],file_list[j][2])
            j+=1
        elif file_list[j][0]=='insu':
            x+=dic[file_list[j][0]]
            ax = fig.add_axes([0.1,1-(x/l)-x/(4*l),0.8,(1/l)*0.65])
            draw_insulation(ax,file_list[j][2],chrs,start,end,file_list[j][1])
            j+=1
        elif file_list[j][0]=='delta_epi':
            x+=dic[file_list[j][0]]
            ax = fig.add_axes([0.1,1-(x/l)-x/(4*l),0.8,(1/l)*1.5])
            y=int(file_list[j][1])-1
            z=int(file_list[j][2])-1
            draw_diff_epigenetic(file_list[y][3],file_list[z][3],ax,chrs,start,end,file_list[j][3],file_list[j][4],file_list[j][5],'none')
            j+=1
        elif file_list[j][0]=='gtf':
            x+=dic[file_list[j][0]]
            ax = fig.add_axes([0.1,1-(x/l)-x/(4*l),0.8,(1/l)*3.0])
            gtf_bed  = filter_gtffile(file_list[j][1],chrs,start,end,args.genetype)
            draw_genes(ax,gtf_bed,start,end)
            #gtf_bed.to_csv(args.output,sep='\t',header=None,index=False)
            j+=1
        elif file_list[j][0]=='DI':
            x+=dic[file_list[j][0]]
            ax = fig.add_axes(fig.add_axes([0.1,1-(x/l)-x/(4*l),0.8,(1/l)*0.9]))
            draw_bar(ax,file_list[j][1],chrs,start,end,1000,-1000)
            j+=1
        elif file_list[j][0]=='AB':
            x+=dic[file_list[j][0]]
            ax = fig.add_axes(fig.add_axes([0.1,1-(x/l)-x/(4*l),0.8,(1/l)*0.7]))
            draw_AB(file_list[j][3],int(file_list[j][1]),chrs,start,end,file_list[j][2],ax)
            j+=1
        elif  file_list[j][0]== "boundary":
            title = file_list[j][1]
            x+=dic[file_list[j][0]]
            ax = fig.add_axes([0.1,1-(x/l)+((1/l)*0.45)*1-x/(4*l),0.8,(1/l)*0.47])
            #ax2 = fig.add_axes([0.1,1-(x/l)+((1/l)*0.45)*3-x/(4*l),0.8,(1/l)*0.47])
            #draw_border(ax1,ax2,file_list[j][2],chrs,start,end)
            draw_compare_border(title,ax,file_list[j][3],chrs,start,end,file_list[j][2])
            j+=1
        elif  file_list[j][0]== "CDB":
            title = file_list[j][1]
            x+=dic[file_list[j][0]]
            ax = fig.add_axes([0.1,1-(x/l)+((1/l)*0.45)*1-x/(4*l),0.8,(1/l)*0.47])
            draw_CDB(title,ax,file_list[j][2],chrs,start,end)
            j+=1
    fig.savefig(args.output)


if __name__=="__main__":
    main()


"""
elif file_list[j][0]=='compartment':
    ax = fig.add_axes([0.1,1-j/n,0.8,(1/n)*0.9])
    draw_compartment(ax,file_list[j][1],file_list[j][2],file_list[j][3],x)
    j+=1
    x+=diff_dic[zfile_list[zj][0]]
        
elif file_list[j][0]=='boundary':
    ax = fig.add_axes([0.1,1-(j/n+(1/(2*n))),0.8,(1/n)/4])
    draw_insu(ax,file_list[j][1],file_list[j][2],file_list[j][3],x)
    j+=1
    x+=diff_dic[file_list[j][0]]
        
"""

