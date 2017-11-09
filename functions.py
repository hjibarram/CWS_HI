#!/usr/bin/python
import sys
import numpy as np
#from numpy import float_
#from numpy import absolute as abs
#from numpy import random as ran
import time
import pyfits as pyf
#from pyfits import getheader as ghead
from pyfits import getdata as gdata
#from pyfits import writeto as wfits
#import os.path as ptt
#from scipy.interpolate.interpolate import interp1d
#from scipy.interpolate import splev, splrep
import os.path as pt
import matplotlib
#import imp
#import cosmolopy.distance as cd
#from pylab import *
import math
#import illustris_python as il
#import h5py
#import random as rat

def get_ages(name,fd=0,dir_data="./"):
    if pt.exists(dir_data+"/"+name+".Pipe3D.cube.fits.gz") == False or fd == 1:
        get_data(name,path_bin=path_bin,dir_o=dir_data,fd=fd,ver=ver)
    file=dir_data+"/"+name+".Pipe3D.cube.fits.gz"
    [pdl_cube, hdr]=gdata(file,2, header=True)
    ages_l,met_l,ages2_l=get_ages_met(hdr)
    ages=np.log10(ages_l)+9
    n_age=len(ages)
    d_ages=np.zeros(n_age)
    for i in range(0, n_age):
        if i == 0:
            age_s=10.0**((ages[i]+ages[i+1])/2.0)
            age_i=0.0
        elif i == n_age-1:
            age_i=10.0**((ages[i]+ages[i-1])/2.0)
            age_s=2.0*10.0**(ages[i])-age_i
        else:
            age_i=10.0**((ages[i]+ages[i-1])/2.0)
            age_s=10.0**((ages[i]+ages[i+1])/2.0)
        d_ages[i]=np.abs(age_s-age_i)
    return ages_l,d_ages

def def_index_sfh(it1=-1,it2=-1,n_met=4,n_age=39):
    na=n_age*n_met
    nb=n_age*(n_met+1)
    #print it1,it2
    if it1 == -1 and it2 == -1:
        nf=-1
    if it1 > -1 and it2 == -1:
        nf=na+it1
    if it1 == -1 and it2 > -1:
        nf=nb+it2
    if it1 > -1 and it2 > -1:
        nf=it1*n_met+it2
        if nf > na:
            nf=-1
    return nf[0]

def asosiate_val_s(val,val_s):
    val_a=[]
    ns=len(val_s)
    ban=0
    val_t=sorted(val_s)
    val_a.extend([val_t[0]])
    ind_s=np.zeros(1,dtype=np.int)
    ind_s[:]=-100
    for i in range(1, ns):
        if val_t[i-1] > val_t[i]:
            ban =1
        if val_t[i-1] < val_t[i] and ban == 0:
            val_a.extend([val_t[i]])
    n_val=len(val_a)
    ind_val=val_definition_l(val,val_a)
    for i in range(0, n_val):
        if len(ind_val[i]) > 0:
            nt=np.where((val_s == val_a[i]))[0]
            ind_s[ind_val[i]]=nt[0]
    return ind_s

def val_definition_l(val,val_ssp):
    val_l=val
    n_val=len(val_ssp)
    ind=[]
    for i in range(0, n_val):
        if i < n_val-1:
            dval=(val_ssp[i+1]+val_ssp[i])/2.
        else:
            dval=val_ssp[i]+1.0
        if i == 0:
            val1=0
        else:
            val1=val2
        val2=dval
        nt=np.where((val_l >= val1) & (val_l <= val2))
        ind.extend([nt[0]])
    return ind

def get_ages_met(hdr,n_met=4,n_age=39):
    nt=n_age*(n_met+1)+n_met
    n_t=n_age*n_met
    a_age=[]
    a_agel=[]
    a_met=[]
    for i in range(0, n_t):
        if (i % n_met) == 0:
            a_agel.extend([np.float(hdr["DESC_"+str(i)].replace("Luminosity Fraction for age-met ","").replace(" SSP","").split("-")[0])])
    for i in range(n_t, nt):
        if i < n_age*(n_met+1):
            a_age.extend([np.float(hdr["DESC_"+str(i)].replace("Luminosity Fraction for age ","").replace(" SSP",""))])
        else:
            a_met.extend([np.float(hdr["DESC_"+str(i)].replace("Luminosity Fraction for met ","").replace(" SSP",""))])
    a_age=np.array(a_age)
    a_agel=np.array(a_agel)
    a_met=np.array(a_met)
    return a_age,a_met,a_agel

def limits(yt,low_v=-1.0,up_v=1.0,logt=0):
    import warnings
    warnings.filterwarnings("ignore")
    yt[np.where(yt == 0)]=np.nan
    upvalue=math.ceil(np.nanmax(yt)/.05)*.05
    #print upvalue
    if np.isinf(upvalue):
        upvalue=up_v
    #lovalue=math.ceil(np.nanmin(yt)/.05)*.05
    lovalue=(np.nanmin(yt)/.05)*.05
    if np.isinf(lovalue):
        lovalue=low_v
    #if lovalue
    if logt == 1:
        lovalue=np.log10(lovalue)
        upvalue=np.log10(upvalue)
        if np.isinf(lovalue):
            lovalue=upvalue-1.5
    #if upvalue*.5 <= lovalue:
    #    if upvalue*.5 > lovalue:
    #        lovalue=upvalue-1
    #else:
    #    lovalue=low_v
    #print lovalue,upvalue
    return lovalue,upvalue

def wfits(name, data, hdr):
    if pt.exists(name) == False:
        wfit(name,data,hdr)
    else:
        name1=name.replace("\ "," ")
        name1=name1.replace(" ","\ ")
        sycall("rm "+name1)
        wfit(name,data,hdr)

def sycall(comand):
    import os
    linp=comand
    os.system(comand)
    
def get_data(name,ver="2.1.2",dir_o="./",path_bin="",fd=0):
    vr="v"+ver.replace(".","_")
    nam=name.split('-')
    dir="https://data.sdss.org/sas/dr14/manga/spectro/pipe3d/"+vr+"/"+ver+"/"+nam[1]+"/"+name+".Pipe3D.cube.fits.gz"
    if pt.exists(dir_o) == False:
        dir_t=dir_o
        DIRS=dir_t.split("/")
        DRT=""
        for DR in DIRS:
            DRT=DRT+DR+"/"
            call="mkdir -p "+DRT
            sycall(call)
    call=path_bin+"wget "+dir
    if pt.exists(dir_o+"/"+name+".Pipe3D.cube.fits.gz") == False or fd == 1:
        sycall(call)
        time.sleep(2)
        if pt.exists(name+".Pipe3D.cube.fits.gz") == False:
            print "There is no file name"
        else:
            print "Download successfully"
            if dir_o is not "./":
                call="mv "+name+".Pipe3D.cube.fits.gz "+dir_o+"/"+name+".Pipe3D.cube.fits.gz"
                sycall(call)
    else:
        print "The file already exists"
    
def map_plot(image,ages_l,overc=np.array([-1]),color_b=0,xli=0,yli=0,legen='',pdf=0,form='pdf',dir='',fname='map',title='map',minval=np.nan,maxval=np.nan,logt=0,pix=0.5,blackw=0,rgb=0,rad=1,no_cb=0):
    import warnings
    warnings.filterwarnings("ignore")
    image[np.isnan(image)]=0
    #print np.isnan(minval)
    #print np.isnan(maxval)
    if np.isnan(minval) == True or np.isnan(maxval) == True:
        #print "A"
        minval,maxval=limits(image,logt=logt)
    if logt == 1:
        mapf=np.log10(image)
        textl1='log['
        textl2=']'
    else:
        mapf=image
        textl1=''
        textl2=''
    if overc.shape[0] == 1:
        over=0
    else:
        over=1
    nx=image.shape[0]
    ny=image.shape[1]
    ft=open("ctable", "r")
    nc=256
    g1=np.zeros(nc)
    r1=np.zeros(nc)
    b1=np.zeros(nc)
    l=np.zeros(nc)
    for cmap in ft:
        data=cmap.split(" ")
        data=filter(None,data)
        nc=int(float(data[0]))
        r1[nc-1]=(float(data[1]))/255.
        g1[nc-1]=(float(data[2]))/255.
        b1[nc-1]=(float(data[3]))/255.
        l[nc-1]=nc/255.
    bright=1 
    contrast=0.5
    r1[0]=1.0
    g1[0]=1.0
    b1[0]=1.0        
    nc=256
    my_rgb=np.zeros([nc,3])
    my_rgb[:,0]=r1
    my_rgb[:,1]=g1
    my_rgb[:,2]=b1
    x1=np.arange(0,nx,1)-nx/2+1
    y1=np.arange(0,ny,1)-ny/2+1
    [X,Y]=np.meshgrid(x1*pix,y1*pix)

    my_cmap = matplotlib.colors.ListedColormap(my_rgb, name='my_name')
    #if pdf == 1:
    #    matplotlib.use('Agg')
    from matplotlib import cm
    if "vel" in fname:
        my_cmap=cm.rainbow
    if "age" in fname:
        my_cmap=cm.rainbow
    if blackw == 1:
        my_cmap=cm.get_cmap('Greys')
    if "V" in fname:
        my_cmap = matplotlib.colors.ListedColormap(my_rgb, name='my_name')      
    if color_b == 1:
        my_cmap = matplotlib.colors.ListedColormap(my_rgb, name='my_name')
    my_cmap.set_under("white")
    if rgb == 1:
        my_cmap=None
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    font0 = FontProperties()
    fig= plt.figure(figsize=(6,5.5))
    ax = fig.add_axes([0.08, 0.11, 0.78, 0.82])
    plt.xlabel(r"$\alpha\ [arcsec]$",fontsize=14)
    plt.ylabel(r"$\delta\ [arcsec]$",fontsize=14)
    ax.yaxis.set_label_coords(-0.05, 0.5)
    plt.title(title,fontsize=15)
    font = font0.copy()
    font.set_weight('bold')
    if xli > 0:
        ax.set_xlim(-xli,xli)
    if yli > 0:
        ax.set_ylim(-yli,yli)
    plt.text(-nx*pix/2.*0.95,-ny*pix/2*0.95,ages_l, fontsize=15)
    font.set_weight('normal')
    if rad == 1:
        levels=[0.5,1,1.5]
        labels = ['$%3.1f' % val + " R_{50}%$" for val in levels ]
        fmt={}
        for l,s in zip( levels, labels ):
            fmt[l] = s
    else:
        levels=15
        if over == 1:
            overc[np.where(np.isfinite(overc) == False)]=0
            max_l=np.amax(overc[np.where(overc > 0)])
            min_l=np.amin(overc[np.where(overc > 0)])
            min_l=min_l+(max_l-min_l)*0.03
            levels=np.logspace(np.log10(min_l),np.log10(max_l),6)

    im=ax.imshow(mapf,cmap=my_cmap,vmin=minval, vmax=maxval, interpolation='nearest', origin='lower', extent=[-nx*pix/2,nx*pix/2,-ny*pix/2,ny*pix/2])
    if over == 1:
        CS = ax.contour(overc, levels, linewidths=2, extent=(-nx*pix/2,nx*pix/2,-ny*pix/2,ny*pix/2), colors=('black','black','black'))
        if rad == 1:
            ax.clabel(CS, levels, inline=1,fmt=fmt,fontsize=14)
    if rgb != 1 or no_cb !=1:
        cbar_ax = fig.add_axes([0.85, 0.11, 0.03, 0.82])
        cbf=fig.colorbar(im, cax=cbar_ax)
        #cbf.set_label('$'+textl1+legen+textl2+'$',fontsize=15)
        cbf.set_label(textl1+legen+textl2,fontsize=15)
    if pdf == 1:
        plt.savefig(dir+fname+'.'+form)
    else:
        plt.show()
    plt.close()
    
def plot_maps(name,type="Ha_flux",dir_p="./",hdu=-1,channel=-1,dir_data="./",path_bin="",ver="2.1.2",fd=0,plot=1,outp=0,pixl=0.5,blackw=0,log_f=0):
    import warnings
    warnings.filterwarnings("ignore")
    if pt.exists(dir_p) == False:
        dir_t=dir_p
        DIRS=dir_t.split("/")
        DRT=""
        for DR in DIRS:
            DRT=DRT+DR+"/"
            call="mkdir -p "+DRT
            sycall(call)
    ban=0
    bta=0
    if "help" in type or "Help" in type or "HELP" in type:
        ban=3
        print "Usage of plot_maps function:"
        print "plot_maps(manga_name,type=TYPE,dir_p=dip,hdu=HDU,channel=CHANEL,dir_data=dir_d,plot=plot,outp=outp,log_f=log_f)"
        print ""
        print "type: name of the map that you want to plot, default=Ha_flux"
        print "hdu: number of the extension that you want to use, if type is given hdu is not mandatory"
        print "channel: number of the channel that you want to use, if type is given channel is not mandatory"
        print "dir_p: Directory where you want to save the output plots, default is ./"
        print "dir_p: Directory where you have the input fit files, default is ./, if there is no files it will be automatically download"
        print "plot: If 0 there are no output plots, if 1 there are a pdf output plots, if 2 there are screen output plots, default 1"
        print "outp: If 1, the numerical array of the map TYPE will be given as an output, default 0"
        print "log_f: If 0 the flux plots will be plotted in lineal form, If 1 the flux plots will be plotted in logarithmic form, default 0"
        print ""
        print "These are the possible values of hdu, channel and type variables:"
        print ""
        print "HDU,CHANNEL,TYPE comments"
        print "1,0,V-band"
        print "1,1,Continuum_seg_SSP"
        print "1,2,Continuum_dez"
        print "1,3,Median_flux"
        print "1,4,StdDev_flux"
        print "1,5,Age_LW"
        print "1,6,Age_MW"
        print "1,7,Error_age"
        print "1,8,Metallicity_LW (where Z=0.02 is solar metallicity)"
        print "1,9,Metallicity_MW"
        print "1,10,Error_met"
        print "1,11,Dust_attenuation (Av in mag)"
        print "1,12,Error_dust_attnuation"
        print "1,13,Velocity_map_stellar"
        print "1,14,Error_velocity_stellar"
        print "1,15,Velocity_dispersion_stellar"
        print "1,16,Error_velocity_dispersion_stellar"
        print "1,17,Mass-to-light_ratio_stellar (log10)"
        print "1,18,Stellar_Mass  (log10(Msun)/spaxels^2 without dust correction)"
        print "1,19,Stellar_Mass_d (log10(Msun)/spaxels^2 with dust correction)"
        print ""
        print "For the analysis of the Stellar decomposition the type values are"
        print "HDU,CHANNEL,TYPE comments"
        print "2,X,SFH_AGE_MET (or SFH1_AGE or SFH2_MET)"
        print "The values of AGE and MET (or only AGE or only MET) depends on the desired map that the user wants"
        print "AGE must be given in Gyr and MET (where Z=0.02 is solar metallicity)"
        print "The values of AGE and MET are not necessary to be equal to the SSP library, the code automatically selects the nearest AGE"
        print "or MET of the SSP library"
        print  "EXAMPLE: SFH_14_"
        print "The values of X range from 0 to 198, see https://data.sdss.org/datamodel/files/MANGA_PIPE3D/MANGADRP_VER/PIPE3D_VER/PLATE/manga.Pipe3D.cube.html for more info."
        print ""
        print "For the emision lines the type values are Y_TYPE, the channel values are i+X, with Y and X as:" 
        print "Flux:  (X=0)    to get the flux intensity of the i analyzed emission line"
        print "Velm:  (X=57)   to get the velocity map of the i analyzed emission line" 
        print "Disp:  (X=114)  to get the velocity dispersion map of the i analyzed emission line"
        print "Ew:    (X=171)  to get the equivalent width of the i analyzed emission line"   
        print "Eflux: (X=228)  to get the error estimated for the flux intensity of the i analyzed emission line"
        print "Evel:  (X=285)  to get the error estimated for the velocity map of the i analyzed emission line"
        print "Edisp: (X=342)  to get the error estimated for the velocity dispersion of the i analyzed emission line"
        print "Eew:   (X=399)  to get the eError estimated for the equivalent width of the i analyzed emission line"   
        print "EXAMPLE: The map of Ha equivalent width has a TYPE value of Ew_Ha and a CHANNE value of 45+171=216"
        print "The values of the emission lines are:"
        print ""
        print "HDU,CHANNEL,TYPE comments"
        print "3,0+X, Y_[OII]3727   (3727.4 A)"
        print "3,1+X, Y_H12         (3750 A)"  
        print "3,2+X, Y_H11         (3771 A)"
        print "3,3+X, Y_H10         (3798 A)" 
        print "3,4+X, Y_HeI3819     (3819.40 A)"  
        print "3,5+X, Y_H9          (3835 A)"
        print "3,6+X, Y_[NeIII]3869 (3869 A)"   
        print "3,7+X, Y_H8          (3889 A)"
        print "3,8+X, Y_[NeIII]3967 (3967 A)" 
        print "3,9+X, Y_He          (3970.1 A)"
        print "3,10+X,Y_HeI4026     (4026.29 A)"    
        print "3,11+X,Y_[SII]4069   (4069.17 A)"  
        print "3,12+X,Y_[SII]4077   (4076.72 A)"
        print "3,13+X,Y_Hd          (4101.74 A)"
        print "3,14+X,Y_[FeII]4277  (4276.83 A)"
        print "3,15+X,Y_[FeII]4287  (4287.40 A)"
        print "3,16+X,Y_[FeII]4320  (4319.62 A)"
        print "3,17+X,Y_Hg          (4340.47 A)"
        print "3,18+X,Y_[OIII]4363  (4363.21 A)"
        print "3,19+X,Y_[FeII]4414  (4413.78 A)"
        print "3,20+X,Y_[FeII]4416  (4416.27 A)"
        print "3,21+X,Y_HeI4471     (4471 A)"
        print "3,22+X,Y_[FeIII]     (4657.93 A)"    
        print "3,23+X,Y_HeII        (4686.0 A)"
        print "3,24+X,Y_HeI4713     (4712.98 A)"   
        print "3,25+X,Y_HeI4922     (4922.16 A)"
        print "3,26+X,Y_[OIII]5007  (5006.84 A)"  
        print "3,27+X,Y_[OIII]4959  (4958.91 A)"
        print "3,28+X,Y_Hb          (4861.32 A)"
        print "3,29+X,Y_[FeII]4890  (4889.62 A)"  
        print "3,30+X,Y_[FeII]4905  (4905.34 A)"
        print "3,31+X,Y_[FeII]5112  (5111.6299 A)"
        print "3,32+X,Y_[FeII]5159  (5158.7798 A)"
        print "3,33+X,Y_[NI]        (5199.6001 A)"
        print "3,34+X,Y_[FeII]5262  (5261.6201 A)"
        print "3,35+X,Y_[ClIII]5518 (5517.71 A)"
        print "3,36+X,Y_[ClIII]5538 (5537.60 A)"
        print "3,37+X,Y_OI          (5554.94 A)"
        print "3,38+X,Y_[OI]5577    (5577.3101 A)"
        print "3,39+X,Y_[NII]5755   (5754.52 A)"
        print "3,40+X,Y_HeI5876     (5875.62 A)"
        print "3,41+X,Y_[OI]6300    (6300.2998 A)"
        print "3,42+X,Y_[SIII]6312  (6312.40 A)"
        print "3,43+X,Y_SiII        (6347.28 A)"
        print "3,44+X,Y_[OI]6364    (6363.7798 A)"
        print "3,45+X,Y_Ha          (6562.68 A)"
        print "3,46+X,Y_[NII]6584   (6583.41 A)"
        print "3,47+X,Y_[NII]6548   (6548.08 A)"
        print "3,48+X,Y_HeI6678     (6677.97 A)"
        print "3,49+X,Y_[SII]6717   (6716.39 A)"
        print "3,50+X,Y_[SII]6731   (6730.74 A)"
        print "3,51+X,Y_[ArIII]7136 (7136 A)"
        print "3,52+X,Y_[OII]7325   (7325 A)"
        print "3,53+X,Y_[ArIII]7751 (7751 A)"
        print "3,54+X,Y_[SIII]9067  (9068.6 A)"
        print "3,55+X,Y_[SIII]9531  (9530.6 A)"
        print ""
        print "For the index maps the type values are:"
        print "HDU,CHANNEL,TYPE comments (Blue Range  Red Range Index Range)"
        print "4,0, Hd_index        (4083.500-4122.250    4041.600-4079.750    4128.500-4161.000)"
        print "4,1, Hb_index        (4847.875-4876.625    4827.875-4847.875    4876.625-4891.625)"
        print "4,2, Mgb index       (5160.125-5192.625    5142.625-5161.375    5191.375-5206.375)"
        print "4,3, Fe5270_index    (5245.650-5285.650    5233.150-5248.150    5285.650-5318.150)"
        print "4,4, Fe5335_index    (5312.125-5352.125    5304.625-5315.875    5353.375-5363.375)"
        print "4,5, D4000           (4050.000-4250.000    3750.000-3950.000)"
        print "4,6, Hdmod_index     (4083.500-4122.250    4079.000-4083.000    4128.500-4161.000)"
        print "4,7, Hg_index        (4319.750-4363.500    4283.500-4319.750    4367.250-4419.750)"
        print "4,8, EHd_index       (Error Values)"
        print "4,9, EHb_index       (Error Values)"
        print "4,10,EMgb_index      (Error Values)"
        print "4,11,EFe5270_index   (Error Values)"
        print "4,12,EFe5335_index   (Error Values)"
        print "4,13,ED4000          (Error Values)"
        print "4,14,EHdmod_index    (Error Values)"
        print "4,15,EHg_index       (Error Values)"
        print "Check https://data.sdss.org/datamodel/files/MANGA_PIPE3D/MANGADRP_VER/PIPE3D_VER/PLATE/manga.Pipe3D.cube.html for more info."
    
    if channel >= 0 and hdu == -1:
        print "Please input the HDU value, set type=Help for more info"
    if channel == -1 and hdu >= 0:
        print "Please input the CHANEL value, set type=Help for more info"
    if channel >= 0 and hdu >= 0:
        channel_f=channel
        hdu_f=hdu
        if hdu == 1:
            if channel > 19:
                print "Please input a CHANEL valid value, set type=Help for more info"
            else:
                ban=1
        if hdu == 2:
            if channel > 198:
                print "Please input a CHANEL valid value, set type=Help for more info"
            else:
                ban=1
                bta=1
                hdu_f=hdu
                #print hdu,channel
        if hdu == 3:
            if channel > 454:
                print "Please input a CHANEL valid value, set type=Help for more info"
            else:
                ban=1
        if hdu == 4:
            if channel > 15:
                print "Please input a CHANEL valid value, set type=Help for more info"
            else:
                ban=1
        if hdu > 4 or hdu == 0:
            print "Please input a HDU valid value, set type=Help for more info"
    else:
        if hdu < -1 or channel < -1:
            if hdu < -1:
                print "Please input a HDU valid value, set type=Help for more info"
            if channel < -1:
                print "Please input a CHANNEL valid value, set type=Help for more info"
        else:
            if "V-band" in type:
                channel_f=0
                hdu_f=1
                ban=1 
            if "Continuum_seg_SSP" in type:
                channel_f=1
                hdu_f=1
                ban=1 
            if "Continuum_dez" in type:
                channel_f=2
                hdu_f=1
                ban=1 
            if "Median_flux" in type:
                channel_f=3
                hdu_f=1
                ban=1
            if "StdDev_flux" in type:
                channel_f=4
                hdu_f=1
                ban=1 
            if "Age_LW" in type:
                channel_f=5
                hdu_f=1
                ban=1
            if "Age_MW" in type:
                channel_f=6
                hdu_f=1
                ban=1 
            if "Error_age" in type:
                channel_f=7
                hdu_f=1
                ban=1
            if "Metallicity_LW" in type:
                channel_f=8
                hdu_f=1
                ban=1 
            if "Metallicity_MW" in type:
                channel_f=9
                hdu_f=1
                ban=1
            if "Error_met" in type:
                channel_f=10
                hdu_f=1
                ban=1 
            if "Dust_attenuation" in type:
                channel_f=11
                hdu_f=1
                ban=1
            if "Error_dust_attnuation" in type:
                channel_f=12
                hdu_f=1
                ban=1 
            if "Velocity_map_stellar" in type:
                channel_f=13
                hdu_f=1
                ban=1      
            if "Error_velocity_stellar" in type:
                channel_f=14
                hdu_f=1
                ban=1 
            if "Velocity_dispersion_stellar" in type:
                channel_f=15
                hdu_f=1
                ban=1
            if "Error_velocity_dispersion_stellar" in type:
                channel_f=16
                hdu_f=1
                ban=1 
            if "Mass-to-light_ratio_stellar" in type:
                channel_f=17
                hdu_f=1
                ban=1
            if "Stellar_Mass" in type:
                channel_f=18
                hdu_f=1
                ban=1 
            if "Stellar_Mass_d" in type:
                channel_f=19
                hdu_f=1
                ban=1
            #print type
            bt=0
            if "Flux" in type:
                Xf=0
                bt=1
            if "Velm" in type:
                Xf=57
                bt=1
            if "Disp" in type:
                Xf=114
                bt=1
            if "Ew" in type:
                Xf=171
                bt=1
            if "Eflux" in type:
                Xf=228
                bt=1
            if "Evel" in type:
                Xf=285
                bt=1
            if "Edisp" in type:
                Xf=342
                bt=1
            if "Eew" in type:
                Xf=399
                bt=1
            if "[OII]3727" in type and bt == 1:
                channel_f=0+Xf
                hdu_f=3
                ban=1
            if "H12" in type and bt == 1:
                channel_f=1+Xf
                hdu_f=3
                ban=1
            if "H11" in type and bt == 1:
                channel_f=2+Xf
                hdu_f=3
                ban=1
            if "H10" in type and bt == 1:
                channel_f=3+Xf
                hdu_f=3
                ban=1
            if "HeI3819" in type and bt == 1:
                channel_f=4+Xf
                hdu_f=3
                ban=1
            if "H9" in type and bt == 1:
                channel_f=5+Xf
                hdu_f=3
                ban=1
            if "[NeIII]3869" in type and bt == 1:
                channel_f=6+Xf
                hdu_f=3
                ban=1  
            if "H8" in type and bt == 1:
                channel_f=7+Xf
                hdu_f=3
                ban=1
            if "[NeIII]3967" in type and bt == 1:
                channel_f=8+Xf
                hdu_f=3
                ban=1 
            if "He" in type and bt == 1:
                channel_f=9+Xf
                hdu_f=3
                ban=1
            if "HeI4026" in type and bt == 1:
                channel_f=10+Xf
                hdu_f=3
                ban=1
            if "[SII]4069" in type and bt == 1:
                channel_f=11+Xf
                hdu_f=3
                ban=1
            if "[SII]4077" in type and bt == 1:
                channel_f=12+Xf
                hdu_f=3
                ban=1
            if "Hd" in type and bt == 1:
                channel_f=13+Xf
                hdu_f=3
                ban=1
            if "[FeII]4277" in type and bt == 1:
                channel_f=14+Xf
                hdu_f=3
                ban=1
            if "[FeII]4287" in type and bt == 1:
                channel_f=15+Xf
                hdu_f=3
                ban=1
            if "[FeII]4320" in type and bt == 1:
                channel_f=16+Xf
                hdu_f=3
                ban=1
            if "Hg" in type and bt == 1:
                channel_f=17+Xf
                hdu_f=3
                ban=1
            if "[OIII]4363" in type and bt == 1:
                channel_f=18+Xf
                hdu_f=3
                ban=1
            if "[FeII]4414" in type and bt == 1:
                channel_f=19+Xf
                hdu_f=3
                ban=1
            if "[FeII]4416" in type and bt == 1:
                channel_f=20+Xf
                hdu_f=3
                ban=1
            if "HeI4471" in type and bt == 1:
                channel_f=21+Xf
                hdu_f=3
                ban=1
            if "[FeIII]" in type and bt == 1:
                channel_f=22+Xf
                hdu_f=3
                ban=1
            if "HeII" in type and bt == 1:
                channel_f=23+Xf
                hdu_f=3
                ban=1
            if "HeI4713" in type and bt == 1:
                channel_f=24+Xf
                hdu_f=3
                ban=1
            if "HeI4922" in type and bt == 1:
                channel_f=25+Xf
                hdu_f=3
                ban=1
            if "[OIII]5007" in type and bt == 1:
                channel_f=26+Xf
                hdu_f=3
                ban=1
            if "[OIII]4959" in type and bt == 1:
                channel_f=27+Xf
                hdu_f=3
                ban=1
            if "Hb" in type and bt == 1:
                channel_f=28+Xf
                hdu_f=3
                ban=1
            if "[FeII]4890" in type and bt == 1:
                channel_f=29+Xf
                hdu_f=3
                ban=1
            if "[FeII]4905" in type and bt == 1:
                channel_f=30+Xf
                hdu_f=3
                ban=1
            if "[FeII]5112" in type and bt == 1:
                channel_f=31+Xf
                hdu_f=3
                ban=1
            if "[FeII]5159" in type and bt == 1:
                channel_f=32+Xf
                hdu_f=3
                ban=1
            if "[NI]" in type and bt == 1:
                channel_f=33+Xf
                hdu_f=3
                ban=1
            if "[FeII]5262" in type and bt == 1:
                channel_f=34+Xf
                hdu_f=3
                ban=1
            if "[ClIII]5518" in type and bt == 1:
                channel_f=35+Xf
                hdu_f=3
                ban=1
            if "[ClIII]5538" in type and bt == 1:
                channel_f=36+Xf
                hdu_f=3
                ban=1
            if "OI" in type and bt == 1:
                channel_f=37+Xf
                hdu_f=3
                ban=1
            if "[OI]5577" in type and bt == 1:
                channel_f=38+Xf
                hdu_f=3
                ban=1
            if "[NII]5755" in type and bt == 1:
                channel_f=39+Xf
                hdu_f=3
                ban=1
            if "HeI5876" in type and bt == 1:
                channel_f=40+Xf
                hdu_f=3
                ban=1
            if "[OI]6300" in type and bt == 1:
                channel_f=41+Xf
                hdu_f=3
                ban=1
            if "[SIII]6312" in type and bt == 1:
                channel_f=42+Xf
                hdu_f=3
                ban=1
            if "SiII" in type and bt == 1:
                channel_f=43+Xf
                hdu_f=3
                ban=1
            if "[OI]6364" in type and bt == 1:
                channel_f=44+Xf
                hdu_f=3
                ban=1
            if "Ha" in type and bt == 1:
                channel_f=45+Xf
                hdu_f=3
                ban=1
            if "[NII]6584" in type and bt == 1:
                channel_f=46+Xf
                hdu_f=3
                ban=1
            if "[NII]6548" in type and bt == 1:
                channel_f=47+Xf
                hdu_f=3
                ban=1
            if "HeI6678" in type and bt == 1:
                channel_f=48+Xf
                hdu_f=3
                ban=1
            if "[SII]6717" in type and bt == 1:
                channel_f=49+Xf
                hdu_f=3
                ban=1
            if "[SII]6731" in type and bt == 1:
                channel_f=50+Xf
                hdu_f=3
                ban=1
            if "[ArIII]7136" in type and bt == 1:
                channel_f=51+Xf
                hdu_f=3
                ban=1
            if "[OII]7325" in type and bt == 1:
                channel_f=52+Xf
                hdu_f=3
                ban=1
            if "[ArIII]7751" in type and bt == 1:
                channel_f=53+Xf
                hdu_f=3
                ban=1
            if "[SIII]9067" in type and bt == 1:
                channel_f=54+Xf
                hdu_f=3
                ban=1
            if "[SIII]9531" in type and bt == 1:
                channel_f=55+Xf
                hdu_f=3
                ban=1
            if "Hd_index" in type:
                channel_f=0
                hdu_f=4
                ban=1
            if "Hb_index" in type:
                channel_f=1
                hdu_f=4
                ban=1
            if "Mgb index" in type:
                channel_f=2
                hdu_f=4
                ban=1
            if "Fe5270_index" in type:
                channel_f=3
                hdu_f=4
                ban=1
            if "Fe5335_index" in type:
                channel_f=4
                hdu_f=4
                ban=1
            if "D4000" in type:
                channel_f=5
                hdu_f=4
                ban=1
            if "Hdmod_index" in type:
                channel_f=6
                hdu_f=4
                ban=1
            if "Hg_index" in type:
                channel_f=7
                hdu_f=4
                ban=1
            if "EHd_index" in type:
                channel_f=8
                hdu_f=4
                ban=1
            if "EHb_index" in type:
                channel_f=9
                hdu_f=4
                ban=1
            if "EMgb_index" in type:
                channel_f=10
                hdu_f=4
                ban=1
            if "EFe5270_index" in type:
                channel_f=11
                hdu_f=4
                ban=1
            if "EFe5335_index" in type:
                channel_f=12
                hdu_f=4
                ban=1
            if "ED4000" in type:
                channel_f=13
                hdu_f=4
                ban=1
            if "EHdmod_index" in type:
                channel_f=14
                hdu_f=4
                ban=1
            if "EHg_index" in type:
                channel_f=15
                hdu_f=4
                ban=1
            if "SFH" in type:
                hdu_f=2
                ban=1
            if ban == 0:
                print "Please input a valid TYPE value, set type=Help for more info"
    if ban == 1:
        if pt.exists(dir_data+"/"+name+".Pipe3D.cube.fits.gz") == False or fd == 1:
            get_data(name,path_bin=path_bin,dir_o=dir_data,fd=fd,ver=ver)
        file=dir_data+"/"+name+".Pipe3D.cube.fits.gz"
        [pdl_cube, hdr]=gdata(file,hdu_f, header=True)
        if hdu_f== 1:
            lab1=hdr["DESC_"+str(channel_f)].replace("within the wavelength range","").replace("of the stellar population","")
            com=""
            log_ff=0
            if "flux" in hdr["TYPE_"+str(channel_f)]:
                com="/s/cm^2"
                if log_f == 1:
                    log_ff=1
            if "velocity" in hdr["TYPE_"+str(channel_f)]:
                com="/s"
            lab2=r''+hdr["TYPE_"+str(channel_f)]+' ['+hdr["UNITS_"+str(channel_f)]+com+']'
        if hdu_f == 2:
            #print "A"
            #print bta
            if bta == 0:
                ages_l,met_l,ages2_l=get_ages_met(hdr)
                data=type.split("_")
                if len(data) == 0 or len(data) == 1:
                    print "Please input a valid TYPE value, set type=Help for more info"
                    sys.exit()
                if len(data) == 2:
                    if "SFHa" in data[0]:
                        try:
                            Age=np.float(data[1])
                        except:
                            print "Please input a valid TYPE value, set type=Help for more info"
                            sys.exit()
                        it1=asosiate_val_s(Age, ages_l)
                        channel_f=def_index_sfh(it1=it1)
                    elif "SFHb" in data[0]:
                        try:
                            Met=np.float(data[1])
                        except:
                            print "Please input a valid TYPE value, set type=Help for more info"
                            sys.exit()
                        it2=asosiate_val_s(Met, met_l)
                        channel_f=def_index_sfh(it2=it2)
                    else:
                        print "Please input a valid TYPE value, set type=Help for more info"
                        sys.exit()
                if len(data) == 3:
                    try:
                        Age=np.float(data[1])
                    except:
                        print "Please input a valid TYPE value, set type=Help for more info"
                        sys.exit()
                    try:
                        Met=np.float(data[2])
                    except:
                        print "Please input a valid TYPE value, set type=Help for more info"
                        sys.exit()
                    it2=asosiate_val_s(Met, met_l)
                    it1=asosiate_val_s(Age, ages2_l)
                    channel_f=def_index_sfh(it1=it1,it2=it2)
            if bta == 1:
                channel_f=channel
            log_ff=0
            if log_f == 1:
                log_ff=1
            lab2=r''+hdr["TYPE_"+str(channel_f)]+' [fraction]'
            type='SFH_'+hdr["DESC_"+str(channel_f)].replace(" SSP","").replace("Luminosity Fraction for age-met ","").replace("Luminosity Fraction for age ","").replace("Luminosity Fraction for met ","").replace("-","_")
            lab1=hdr["DESC_"+str(channel_f)].replace("Luminosity Fraction","Lum fraction").replace(" SSP","")
            #sys.exit()
        if hdu_f == 3:
            log_ff=0
            lab1=type.replace("_"," ")
            if "flux" in hdr["NAME"+str(channel_f)]:
                lab2=r'Flux [10^-16 erg/s/cm^2/A]'
                if log_f == 1:
                    log_ff=1
            if "vel" in hdr["NAME"+str(channel_f)]:
                lab2=r'Velocity [km/s]'
            if "disp" in hdr["NAME"+str(channel_f)]:
                lab2=r'Dispersion [km/s]'
            if "EW" in hdr["NAME"+str(channel_f)]:
                lab2=r'Ew [A]'
        if hdu_f == 4:
            log_ff=0
            lab1=type.replace("_"," ").replace("E","Error of ")
            lab2=r''+type.replace("_","").replace("index","").replace("E","e_")+' [A]'
        pdl_map=pdl_cube[channel_f,:,:]
        if outp == 1 and plot == 0:
            return pdl_map
        else:            
            if plot > 0:
                image=pdl_map
                image[np.isnan(image)]=0
                if hdu_f == 1:
                    cont=image
                if hdu_f == 2:
                    cont=image
                if hdu_f == 3:
                    cont=np.array([-1])
                if hdu_f == 4:
                    cont=image
                #if plot == 2:
                #    print lovalue,upvalue
                map_plot(image,lab1,overc=cont,logt=log_ff,legen=lab2,dir=dir_p+"/",pdf=plot,title=name,fname=name+'_map_'+type,pix=pixl,blackw=blackw,rad=0)
            if outp == 1:
                return pdl_map
            
        