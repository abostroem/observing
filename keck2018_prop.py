import numpy as np

from matplotlib import pyplot as plt
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
import supernova_typeII_catalog

def calc_mag_age(date, visibility = None):
    #Assume SN is 100 days old at fall from plateau unless rise time was captured
    start_date = Time('2018-02-01', format='iso')
    end_date = Time('2018-07-01', format='iso')
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


    filename = '{}_keck_target.txt'.format(date.value.split(' ')[0])
    ofile = open(filename, 'w')
    

    for irow in tbdata:
        ofile.write('#########################\n')
        ofile.write('{}\n'.format(irow['SN']))
        ofile.write('#########################\n')
        if visibility is not None:
            if irow['SN'] in visibility:
                is_visible = True
            else:
                is_visible = False
        else: #If no dictionary, the do everything
            is_visible = True
        if is_visible is True:
            if not hasattr(irow['end_plateau_date'], '_mask'): #Skip masked data
                #Calculate the explosion time
                t_expl = Time(irow['end_plateau_date'], format='isot')-t_plateau
                #Calculate the time on the Ni decay tail for the beginning and ending of the semester
                delta_t_tail = date - Time(irow['end_drop_date'], format='isot')
                #Calculate the magnitude at the beginning and ending of the semester
                mag = irow['end_drop_mag'] + (ni_decay_slope * delta_t_tail.value)
            else:
                t_expl = Time(irow['peak/first_date'])-t_rise
                end_drop_date = t_expl+t_rise+t_plateau+t_drop
                end_drop_mag = irow['peak/first_mag']+mag_plateau+mag_drop
                delta_t_tail = date - end_drop_date
                mag = end_drop_mag + (ni_decay_slope * delta_t_tail.value)
            age = date - t_expl
            if age.value > 100:
                #Write to file
                ofile.write('\t -----------------------------\n')
                ofile.write('\t mag on {} = {:2.2f}\n'.format(date.isot.split('T')[0], mag))
                ofile.write('\t age on {} = {:2.2f}\n'.format(date.isot.split('T')[0], age.value))
                ofile.write('\t -----------------------------\n')

    ofile.close()

def check_visibility(date):
    cycler_ls = cycler(color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 
                              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
                              '#bcbd22', '#17becf', '#90EE90', '#808000',
                              'Navy', 'm', '#DAA520', 'k', '#3CB371']*2)
    tbdata = supernova_typeII_catalog.get_cat()
    keck = Observer.at_site('keck')
    visibility = []
    fig = plt.figure(figsize=(11, 8.5))
    ax = fig.add_axes([0.125,0.175,0.6, 0.75])
    ax.axes.set_prop_cycle(cycler_ls)
    night_start = keck.twilight_evening_astronomical(date)
    night_end = keck.twilight_morning_astronomical(date+1*U.day)
    time_arr = Time(np.linspace(night_start.value, night_end.value, 500), format='jd')
    label = []
    for irow in tbdata:
        coords = SkyCoord('{} {}'.format(irow['RA'], irow['DEC']), unit=(U.hourangle, U.deg))
        sn = FixedTarget(name=irow['SN'], coord = coords)
        if np.array(keck.target_is_up(time_arr, sn)).any():
            visibility.append(irow['SN'])
            ax = plot_airmass(sn, keck, time_arr, altitude_yaxis=True, ax=ax)
            label.append(irow['SN'])
    ax.legend(label, bbox_to_anchor=(1.15, 1), fontsize='x-small')
    ax.grid()
    ax.set_title(date.value)
    plt.savefig('{}_keck_visibility.pdf'.format(date.value.split(' ')[0]), orientation='landscape')
    plt.close()
    return visibility

        
        
                
            
if __name__ == "__main__":
    dates = [Time('2017-09-14'), Time('2018-02-01'), Time('2018-07-01')]
    for date in dates:
        visibility = check_visibility(date)
        calc_mag_age(date, visibility = visibility)