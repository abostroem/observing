# Code to help observing planning.

# Example

```
'''
Make an observing plan for Nov 2017 night
'''
import os
import sys
#Clone a local copy from: https://github.com/abostroem/observing
observing_dir = '../../../observing'
if not os.path.exists(observing_dir):
    print('Please download a local copy of the observing scripts from: https://github.com/abostroem/observing')
    sys.exit()
sys.path.append(observing_dir)
import choose_targets
import supernova_typeII_catalog
from matplotlib import pyplot as plt

    
observatory = 'keck'
dates = ['2017-11-17']

#Add objects not in table (e.g. ones that you don't want to calculate the age and current
#magnitude of
sn2016bkv = choose_targets.target('2016bkv', dates[0], ra='10:18:19.31', dec='+41:25:39.3')
iPTF14hls = choose_targets.target('14hls', dates[0], ra='09:20:34.30', dec='+50:41:46.8')

other_objects = [sn2016bkv, iPTF14hls]

for date in dates:
    ofile = open('{}_{}_targets.txt'.format(date, observatory), 'w')
    ofile = choose_targets.write_header(ofile, date, observatory)
    fig = plt.figure(figsize=(11, 8.5))
    ax = fig.add_axes([0.125,0.175,0.6, 0.75])
    ax.axes.set_prop_cycle(choose_targets.cycler_ls)
    label = []
    tbdata = supernova_typeII_catalog.get_cat()
    for snname in tbdata['SN']:
        sn = choose_targets.target(snname, date)
        sn = choose_targets.calc_mag_age(sn)
        sn, ax = choose_targets.check_visibility(sn, observatory, ax)
        choose_targets.write_output_table(ofile, sn)
        if sn.plotted is True:
            label.append(sn.name)
    for sn in other_objects:
        sn.date = date
        sn, ax = choose_targets.check_visibility(sn, observatory, ax)
        if sn.plotted is True:
            label.append(sn.name)
    ofile.close()
    ax = choose_targets.add_airmass_axis(ax)
    ax = choose_targets.add_twilight(ax, observatory, date)
    label.append('12deg twilight')
    ax.legend(label, bbox_to_anchor=(1.10, 1), fontsize='x-small')
    ax.grid()
    ax.set_title('{} {}'.format(date, observatory))
    plt.savefig('{}_{}_visibility.pdf'.format(date, observatory), orientation='landscape')
    plt.close()

```