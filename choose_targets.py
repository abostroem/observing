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

cycler_ls = cycler(color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 
                              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
                              '#bcbd22', '#17becf', '#90EE90', '#808000',
                              'Navy', 'm', '#DAA520', 'k', '#3CB371']*2)

import supernova_typeII_catalog

class target(object):
    def __init__(self,name, date):
        self.name = name
        self.date = Time(date)
        self.mag = None
        self.age = None
        self.visibility = None
        self.visibility_rating = None
        self.plotted = None
        

def calc_airmass(observer, time, target):
    airmass = observer.altaz(time, target).secz
    masked_airmass = np.ma.array(airmass, mask=airmass<1)
    return masked_airmass

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
        if not row['end_plateau_date'].mask: #Skip masked data
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


def check_visibility(sn, observatory, ax):
    tbdata = supernova_typeII_catalog.get_cat()
    observatory = Observer.at_site(observatory)
    night_start = observatory.twilight_evening_astronomical(sn.date)
    night_end = observatory.twilight_morning_astronomical(sn.date+1*U.day)
    time_arr = Time(np.linspace(night_start.value, night_end.value, 500), format='jd')
    indx = tbdata['SN'] == sn.name
    row = tbdata[indx]
    coords = SkyCoord('{} {}'.format(row['RA'].data[0], row['DEC'].data[0]), unit=(U.hourangle, U.deg))
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
            ax = plot_airmass(targ, observatory, time_arr, altitude_yaxis=True, ax=ax)
            sn.plotted = True
        else:
            if (sn.age.value > 75) and (sn.age.value < 550):
                if sn.mag < 24.1:
                    ax = plot_airmass(targ, observatory, time_arr, altitude_yaxis=True, ax=ax)
                    sn.plotted = True
    return sn, ax

def write_output_table(ofile, sn):
    if (sn.age.value > 75) and (sn.age.value < 550):
        if sn.mag < 24.1:
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