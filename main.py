#!/usr/bin/python3

import argparse
import astropy
import matplotlib.pyplot as plt
import numpy as np
import os
import random
import sys
import tweepy
import wikipediaapi as wiki
from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord

# I've added the ability to run the script in different modes. The
# default is to tweet out the results, but you can use the -local
# flag to run things without tweeting, in which case the parameters
# and plots will be displayed. The -manual flag lets you specify
# a particular pulsar by name.
parser = argparse.ArgumentParser()
parser.add_argument('-local', action='store_true')
parser.add_argument('-manual')
args = parser.parse_args()

# Transforms the catalog file into a list that's easier to interpret
with open('{}\\{}'.format(directory, database), 'r') as file:
    catalog = file.read().split('@-----------------------------------------------------------------')[:-1]
    file.close()
for i in range(0, len(catalog)):
    if i == 0:
        catalog[i] = catalog[i].split('\n')[5:][:-1]
    else:
        catalog[i] = catalog[i].split('\n')[1:-1]
    for j in range(len(catalog[i])):
        while '  ' in catalog[i][j]:
            catalog[i][j] = catalog[i][j].replace('  ', ' ')
        if catalog[i][j][-1] == ' ':
            catalog[i][j] = catalog[i][j][:-1]
        catalog[i][j] = catalog[i][j].split()[:2]

# Picks a random pulsar and fetches its key properties
# This block of code could be condensed a bit because some parts
# are repeated between the two sections (manual vs. random)
valid_input = False
while True:
    if valid_input == True:
        break
    # No pulsar has been specified; chooses one at random
    if args.manual == None:
        psr = random.choice(catalog)
        param_names = [psr[i][0] for i in range(len(psr))]
        param_quants = [psr[i][1] for i in range(len(psr))]
        param_dict = {}
        for i in range(len(param_names)):
            param_dict[param_names[i]] = param_quants[i]
        # Some pulsars have a period listed instead of a spin frequency,
        # but most of those have no period derivative, and the vast, vast
        # majority have a spin frequency listed instead.
        psr_name = param_dict['PSRJ']
        if 'PSRB' not in param_dict.keys():
            param_dict['PSRB'] = None
        if '{}.png'.format(psr_name) not in os.listdir('{}\\tweeted_pulsars'.format(directory)) or args.local:
            if set(['PSRJ', 'RAJ', 'DECJ', 'F0', 'F1', 'DM']).issubset(set(param_dict.keys())):
                valid_input = True
    # A pulsar has been specified; checks to see that it's listed
    else:
        for pulsar in catalog:
            param_names = [pulsar[i][0] for i in range(len(pulsar))]
            param_quants = [pulsar[i][1] for i in range(len(pulsar))]
            param_dict = {}
            for i in range(len(param_names)):
                param_dict[param_names[i]] = param_quants[i]
            if 'PSRB' not in param_dict.keys():
                param_dict['PSRB'] = None
            if args.manual in [param_dict['PSRJ'], param_dict['PSRB']]:
                print('{} in catalog.'.format(args.manual))
                valid_input = True
                break
        if valid_input != True:
            print('{} not in catalog. Ending script.'.format(args.manual))
            sys.exit()

psr_bname = param_dict['PSRB']
psr_name = param_dict['PSRJ']
psr_ra = param_dict['RAJ']
psr_dec = param_dict['DECJ']
psr_freq = float(param_dict['F0'])
psr_fdot = float(param_dict['F1'])
psr_period = (1/psr_freq)*u.s
psr_pdot = -(1/psr_freq**2)*psr_fdot
psr_dm = float(param_dict['DM'])*u.pc/(u.cm)**3
psr_char_age = (psr_period/(2*psr_pdot)).to(u.yr)
psr_B_S = 10**12 * np.sqrt(psr_pdot/(10**(-15))) * np.sqrt(psr_period/u.s) * u.G
# It may be worth adding other quantities, like spindown luminosity or distance

# Checks which telescopes should be able to observe the pulsar
# based purely on public declination ranges. It doesn't make the
# distinction between visiblity and detectability, i.e. whether
# a pulsar is too faint to be detected by a given telescope.
visible_telescopes = ''
def dec_to_deg(dec):
    dec = dec.split(':')
    dec = [float(c) for c in dec]
    return Angle(tuple(dec), u.deg).degree
dec_ranges = {'Arecibo': ['-1:00:00', '37:30:00'],
    'CHIME': ['-15:00:00', '90:00:00'],
    'FAST': ['-14:12:00', '65:48:00'],
    'GBT': ['-46:00:00', '90:00:00'],
    'VLA': ['-44:00:00', '90:00:00']}
for telescope in dec_ranges.keys():
    lower = dec_to_deg(dec_ranges[telescope][0])
    upper = dec_to_deg(dec_ranges[telescope][1])
    if lower < dec_to_deg(psr_dec) and dec_to_deg(psr_dec) < upper:
        visible_telescopes += '{}, '.format(telescope)
if visible_telescopes != '':
    visible_telescopes = '\n' + 'Visible from: ' + visible_telescopes
    visible_telescopes = visible_telescopes[:-2]

# Checks to see if Wikipedia page exists under either name and
# adds link to page if one exists
language = "en"
wikipedia = wiki.Wikipedia(language)
pulsar_page = wikipedia.page('PSR_{}'.format(psr_name))
pulsar_b_page = wikipedia.page('PSR_{}'.format(psr_bname))
if pulsar_page.exists():
    wiki_link = '\n' + 'Wikipedia: https://wikipedia.org/wiki/PSR_{}'.format(psr_name)
elif pulsar_b_page.exists():
    wiki_link = '\n' + 'Wikipedia: https://wikipedia.org/wiki/PSR_{}'.format(psr_b_name)
else:
    wiki_link = ''

output = 'Pulsar: {}\n'.format(psr_name) +\
         'RA: {}\n'.format(psr_ra) +\
         'Dec: {}\n'.format(psr_dec) +\
         'Period: {:.3e} s\n'.format(round(psr_period.value, 3)) +\
         'Pdot: {:.3e}\n'.format(psr_pdot) +\
         'DM: {} pc/cm3\n'.format(psr_dm.value) +\
         'Characteristic age: {:.3e} yr\n'.format(psr_char_age.value) +\
         'Surface magnetic field: {:.3e} G'.format(round(psr_B_S.value, 3)) +\
         visible_telescopes +\
         wiki_link

pers = []
pdots = []
names = []
coords = []
# Picks a selection of random pulsars to plot for comparison
for j in range(1000):
    while True:
        psr = random.choice(catalog)
        param_names = [psr[i][0] for i in range(len(psr))]
        param_quants = [psr[i][1] for i in range(len(psr))]
        param_dict = {}
        for i in range(len(param_names)):
            param_dict[param_names[i]] = param_quants[i]
        if set(['PSRJ', 'RAJ', 'DECJ', 'F0', 'F1']).issubset(set(param_dict.keys())):
            name = param_dict['PSRJ']
            if name not in names and name != psr_name:
                names.append(name)
                freq = float(param_dict['F0'])
                fdot = float(param_dict['F1'])
                period = (1/freq)*u.s
                pdot = -(1/freq**2)*fdot
                pers.append(period.value)
                pdots.append(pdot)
                ra_l = param_dict['RAJ'].split(':')
                ra = ra_l[0] + 'h' + ra_l[1] + 'm'
                if len(ra_l) > 2:
                    ra += ra_l[2] + 's'
                dec_l = param_dict['DECJ'].split(':')
                dec = dec_l[0] + 'd' + dec_l[1] + 'm'
                if len(dec_l) > 2:
                    dec += dec_l[2] + 's'
                coords.append([ra, dec])
                break

# P-Pdot diagram
ax1 = plt.subplot(121)
ax1.scatter(pers, pdots)
ax1.scatter([psr_period.value], [psr_pdot], c='r', label=psr_name)
# I currently have Vela, the Crab and Geminga, but I think it might
# be useful to pick well-known pulsars that are more spread out
# across the P-Pdot plane.
comp_psrs = {
    'Vela': [0.0893, 1.250*10**(-13), 'orange'],
    'Crab': [0.0334, 4.204*10**(-13), 'green'],
    'Geminga': [0.2371, 1.097*10**(-13), 'purple']}
for psr in comp_psrs.keys():
    per = comp_psrs[psr][0]
    pdot = comp_psrs[psr][1]
    color = comp_psrs[psr][2]
    ax1.scatter([per], [pdot], c=color, label=psr)

# I wanted to show lines denoting ages as an interesting comparison;
# see Handbook of Pulsar Astronomy, Lorimer & Kramer, Fig. 1.13 for
# what I'm aiming for.
p = np.logspace(-3, 1)
suffix = {5: '00 kyr', 6: ' Myr', 7: '0 Myr', 8: '00 Myr', 9: ' Gyr', 10: '0 Gyr'}
for num in suffix.keys():
    age = 10**num * u.yr
    age = (age.to(u.s)).value
    pdots = p/(2*age)
    ax1.plot(p, pdots, 'k--')
    ax1.text(p[1], 4*pdots[1], '1{}'.format(suffix[num]), rotation=20)

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(10**(-3), 10)
ax1.set_ylim(10**(-22), 10**(-10))
ax1.set_xlabel('Period (s)')
ax1.set_ylabel('Period derivative (s s$^{-1}$)')
ax1.legend()

# Skymap in galactic coordinates. The projection should arguably
# be something besides an Aitoff projection, but eh.
ax2 = plt.subplot(122, projection='aitoff')

l_list = []
b_list = []
for i in range(len(names)):
    c = SkyCoord(ra=coords[i][0], dec=coords[i][1], frame='icrs')
    c = c.galactic
    l_rad = c.l.wrap_at(180*u.deg).radian
    b_rad = c.b.radian
    l_list.append(l_rad)
    b_list.append(b_rad)
ax2.grid(True)
ax2.scatter(l_list, b_list)
    
ra_l = psr_ra.split(':')
ra = ra_l[0] + 'h' + ra_l[1] + 'm'
if len(ra_l) > 2:
    ra += ra_l[2] + 's'
    dec_l = psr_dec.split(':')
    dec = dec_l[0] + 'd' + dec_l[1] + 'm'
    if len(dec_l) > 2:
        dec += dec_l[2] + 's'
psr_c = SkyCoord(ra=ra, dec=dec, frame='icrs')
psr_c = psr_c.galactic
psr_l_rad = psr_c.l.wrap_at(180* u.deg).radian
psr_b_rad = psr_c.b.radian
ax2.scatter(psr_l_rad, psr_b_rad, c='r', label=psr_name)
ax2.legend()

fig = plt.gcf()
fig.set_size_inches(12.0, 7.0)

if args.local:
    print(output)
    plt.show()
else:
    plt.savefig("{}\\tweeted_pulsars\\{}.png".format(directory, psr_name), bbox_inches="tight")
    try:
        auth = tweepy.OAuthHandler(consumer_key, consumer_secret)
        auth.set_access_token(access_token, access_token_secret)
        api = tweepy.API(auth, wait_on_rate_limit=True,
            wait_on_rate_limit_notify=True)
        api.update_with_media(filename="{}\\tweeted_pulsars\\{}.png".format(directory, psr_name), status=output)
        print("Tweeting successful!")
    except tweepy.error.TweepError:
        print("Invalid account credentials.")
