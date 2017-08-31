'''
Make an observing plan for Sept 2017 night and Investigate
possible targets for 2018A semester at Keck
'''

import choose_targets
import supernova_typeII_catalog
from matplotlib import pyplot as plt

    
observatory = 'keck'
dates = ['2017-09-14', '2018-02-01','2018-05-01', '2018-07-01']
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
    ofile.close()
    ax.legend(label, bbox_to_anchor=(1.15, 1), fontsize='x-small')
    ax.grid()
    ax.set_title('{} {}'.format(date, observatory))
    plt.savefig('{}_{}_visibility.pdf'.format(date, observatory), orientation='landscape')
    plt.close()
