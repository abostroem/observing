import numpy as np

from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib import dates

from cycler import cycler

from astropy import units as U
from astropy.time import Time, TimeDelta

from astroplan import Observer
from astropy.coordinates import SkyCoord
from astroplan import FixedTarget
from astroplan.plots import plot_airmass


from supernova import LightCurve


'''
Used for Keck proposal 2018A
'''

cycler_color = cycler(color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 
                              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
                              '#bcbd22', '#17becf', '#90EE90', '#808000',
                              'Navy', 'm', '#DAA520', 'k', '#3CB371'])
cycler_ls = cycler(linestyle=['-', '--', ':'])
combine_cycler = cycler_ls*cycler_color


import supernova_typeII_catalog

class target(object):
    def __init__(self,name, date, age_lim = 550, mag_lim=24, age_min=75., ra=None, dec=None):
        self.name = name
        self.date = Time(date)
        self.mag = None
        self.age = None
        self.visibility = None
        self.visibility_rating = None
        self.plotted = None
        self.age_lim = age_lim
        self.mag_lim = mag_lim
        self.ra = ra
        self.dec = dec
        self.age_min = age_min
        

def calc_airmass(observer, time, target):
    airmass = observer.altaz(time, target).secz
    masked_airmass = np.ma.array(airmass, mask=airmass<1)
    return masked_airmass

def plot_visibility(target, observatory, time, ax=None, fmt={'color':None, 'linestyle':'-'}):
    '''
    time is assumed to be an array of astropy time objects
    observer is an astroplan observer
    target is an astroplan fixed target
    '''
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
    altitude = observatory.altaz(time, target).alt
    ax.plot_date(time.plot_date, altitude, label=target.name, color=fmt['color'], ls=fmt['linestyle'], marker='None')
    
    date_formatter = dates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(date_formatter)
    plt.setp(ax.get_xticklabels(), rotation=30, ha='right', fontsize='x-small')   
    plt.setp(ax.get_yticklabels(), fontsize='x-small') 
    ax.set_ylim(5, 90)
    ax.set_ylabel('Altitude (deg)')
    ax.set_xlabel('Time (UTC)')
    return ax
    
def add_airmass_axis(ax):
    ylim = ax.get_ylim()
    yticks = ax.get_yticks()
    yticks_rad = (yticks*U.deg).to(U.radian).value
    ax_airmass = ax.twinx()
    ax_airmass.set_yticks(yticks)
    ax_airmass.set_ylim(ylim)
    y = 1/np.sin(yticks_rad+0.000001)
    y_txt = ['{:.2f}'.format(itick) for itick in y]
    ax_airmass.set_yticklabels(y_txt, fontsize='x-small') #z = cos(90-alt) = sin(alt)
    ax_airmass.set_ylabel('Airmass')
    ax_airmass.yaxis.set_label_coords(1.10, 0.2)
    ax_airmass.grid(color='WhiteSmoke')  
    return ax

def add_twilight(ax, observatory, date):
    date = Time(date)
    observatory = Observer.at_site(observatory)
    start = observatory.twilight_evening_nautical(date)
    end = observatory.twilight_morning_nautical(date+1*U.day)
    ylim = ax.get_ylim()
    ax.plot_date([start.plot_date, start.plot_date], ylim, 'k:')
    ax.plot_date([end.plot_date, end.plot_date], ylim, 'k:', label='twilight')
    return ax
    

def calc_mag_age(sn):
    #Assume SN is 100 days old at fall from plateau
    t_plateau = 100.*U.day
    ni_decay_slope = .01 #mag/day = 1 mag/100day
    #If no fall from plateau but max is caught then:
        #Assume age is max - 10 days
        #Assume a drop in magnitude of 2 over plateau and 2.5 in fall from plateau then normal decay
        #Assume 100 day plateau
        #Assume 20 day to fall from plateau (a0=1.75, w0=5)
    mag_plateau = 2
    mag_drop = 2.5
    t_rise = 10*U.day
    t_drop = 20*U.day

    tbdata = supernova_typeII_catalog.get_cat()
    indx = tbdata['SN'] == sn.name
    row  = tbdata[indx]
    if sn.visibility is None:
        is_visible = True
    else:
        is_visible = sn.visibility
    if is_visible is True:
        if not hasattr(row['end_plateau_date'],'mask'): #Skip masked data
            #Calculate the explosion time
            t_expl = Time(row['end_plateau_date'], format='isot')-t_plateau
            #Calculate the time on the Ni decay tail for the beginning and ending of the semester
            delta_t_tail = sn.date - Time(row['end_drop_date'], format='isot')
            #Calculate the magnitude at the beginning and ending of the semester
            sn.mag = row['end_drop_mag'] + (ni_decay_slope * delta_t_tail.value)
        else:
            t_expl = Time(row['peak/first_date'])-t_rise
            end_drop_date = t_expl+t_rise+t_plateau+t_drop
            end_drop_mag = row['peak/first_mag']+mag_plateau+mag_drop
            delta_t_tail = sn.date - end_drop_date
            sn.mag = end_drop_mag + (ni_decay_slope * delta_t_tail.value)
        sn.age = sn.date - t_expl
    return sn


def check_visibility(sn, observatory, ax, fmt='-'):
    '''
    RA and DEC should be in this format: 02:35:30.15 -09:21:15.0
    '''
    tbdata = supernova_typeII_catalog.get_cat()
    observatory = Observer.at_site(observatory)
    night_start = observatory.twilight_evening_nautical(sn.date)
    night_end = observatory.twilight_morning_nautical(sn.date+1*U.day)
    time_arr = Time(np.linspace((night_start - 30*U.minute).value, (night_end+30*U.minute).value, 500), format='jd')
    indx = tbdata['SN'] == sn.name
    row = tbdata[indx]
    if sn.ra is None and sn.dec is None:
        coords = SkyCoord('{} {}'.format(row['RA'].data[0], row['DEC'].data[0]), unit=(U.hourangle, U.deg))
    else:
        coords = SkyCoord('{} {}'.format(sn.ra, sn.dec), unit=(U.hourangle, U.deg))
    targ = FixedTarget(name=sn.name, coord = coords)
    if np.array(observatory.target_is_up(time_arr, targ)).any():
        sn.visibility = True
        airmass = calc_airmass(observatory, time_arr, targ)
        if airmass.min() >=3: sn.visibility_rating = 5
        elif airmass.min() >=2: sn.visibility_rating = 4
        elif airmass.min() >=1.5: sn.visibility_rating =3
        elif airmass.min() >1: sn.visibility_rating = 2
        else: sn.visibility_rating = 1
        if (sn.age is None) and (sn.mag is None):
            ax = plot_visibility(targ, observatory, time_arr, ax=ax, fmt=fmt)
            ax.grid(True, color='WhiteSmoke')  
            sn.plotted = True
        else:
            if (sn.age.value > sn.age_min) and (sn.age.value < sn.age_lim):
                if sn.mag < sn.mag_lim:
                    ax = plot_visibility(targ, observatory, time_arr, ax=ax, fmt=fmt)
                    sn.plotted = True    
    return sn, ax

def write_output_table(ofile, sn):
    if (sn.age.value > 75) and (sn.age.value < sn.age_lim):
        if sn.mag < sn.mag_lim:
            if sn.visibility is True:
                ofile.write('#########################\n')
                ofile.write('{}; VISIBILITY={}\n'.format(sn.name, sn.visibility_rating))
                ofile.write('#########################\n')
                ofile.write('\t mag on {} = {:2.2f}\n'.format(sn.date.isot.split('T')[0], sn.mag.data[0]))
                ofile.write('\t age on {} = {:2.2f}\n'.format(sn.date.isot.split('T')[0], sn.age.value[0]))
    return ofile

def write_header(ofile, date, observatory):
    ofile.write('===============================================================\n')
    ofile.write('Targets for the night of {} observred at {}\n\n'.format(date, observatory))
    ofile.write('===============================================================\n')
    ofile.write('\nAssume SN is 100 days old at fall from plateau. End of plateau and end of fall from \n\
plateau are approximated visually from LCO SNEX data. If no fall from plateau \n\
can be measured, assume explosion is max-10 days, a drop of 2 mag over 100 \n\
days in the plateau, and a drop of 2.5 mag over 20 days in the fall from \n\
plateau. Remaining magnitude change is calculated from the nickel decay. \n\
Assume a nickel decay slope of 0.1 mag/day.\n\n')
    ofile.write('Visibility Rating Scale:\n------------------------------\n\t1: airmass =1 \n\t2: 1<airmass<1.5 \n\t3: 1.5<airmass<2\n\t4: 2<airmass<3\n\t5: airmass>3\n\n')
    return ofile