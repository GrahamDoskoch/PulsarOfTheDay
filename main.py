#!/usr/bin/python3

import argparse
import astropy
import matplotlib.pyplot as plt
import numpy as np
import os
import random
import requests
import sys
import tweepy
import wikipediaapi as wiki
from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from bs4 import BeautifulSoup
from urllib.request import urlopen

# The directory where the database is being stored must be specified, as must the relevant credentials
# for the Twitter account. I've removed mine from the public version of the script, for security
# reasons. If you're just running this locally, though, the Twitter information can be ignored.
directory = ''
database = ''
consumer_key = ''
consumer_secret = ''
access_token = ''
access_token_secret = ''

pta_file = 'pta_pulsars.yaml'
pta_database = yaml.load(open(pta_file))

# The script can be run in different modes. The default is to tweet
# out the results, but you can use the -local flag to run things
# without tweeting, in which case the parameters and plots will be
# displayed. The -manual flag lets you specify a particular pulsar
# by name.
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
        param_dict['RM'] = None
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
            param_dict['RM'] = None
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
psr_rm = param_dict['RM']
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
    'Effelsberg': ['-30:00:00', '90:00:00'],
    'FAST': ['-14:12:00', '65:48:00'],
    'GBT': ['-46:00:00', '90:00:00'],
    'Parkes': ['-90:00:00', '25:00:00'],
    'VLA': ['-44:00:00', '90:00:00']}
for telescope in dec_ranges.keys():
    lower = dec_to_deg(dec_ranges[telescope][0])
    upper = dec_to_deg(dec_ranges[telescope][1])
    if lower < dec_to_deg(psr_dec) and dec_to_deg(psr_dec) < upper:
        visible_telescopes += '{}, '.format(telescope)
if visible_telescopes != '':
    visible_telescopes = '\n' + 'Visible from: ' + visible_telescopes
    visible_telescopes = visible_telescopes[:-2]

# Currently a problem with getting the pulsar page
"""
# Checks to see if Wikipedia page exists under either name and
# adds link to page if one exists
language = "en"
wikipedia = wiki.Wikipedia(language)
pulsar_page = wikipedia.page('PSR_{}'.format(psr_name))
pulsar_b_page = wikipedia.page('PSR_{}'.format(psr_bname))
print(pulsar_page)
if pulsar_page.exists():
    print('y')
    wiki_link = '\n' + 'Wikipedia: https://wikipedia.org/wiki/PSR_{}'.format(psr_name)
elif pulsar_b_page.exists():
    wiki_link = '\n' + 'Wikipedia: https://wikipedia.org/wiki/PSR_{}'.format(psr_bname)
else:
    wiki_link = ''
"""
wiki_link = ''

# Looks for flux densities. If values are listed at 400 or 1400 MHz,
# those are used; if they're not but flux densities are listed at other
# frequencies, one of those is randomly chosen.
flux_densities = []
if 'S400' in param_dict.keys():
    flux_freq = 400
    flux_val = param_dict['S400']
    is_flux = True
elif 'S1400' in param_dict.keys():
    flux_freq = 1400
    flux_val = param_dict['S1400']
    is_flux = True
else:
    for key in param_dict.keys():
        if key[0] == 'S' and key[1:].isnumeric():
            flux_densities.append([int(key[1:]), param_dict[key]])
    if flux_densities != []:
        pair = random.choice(flux_densities)
        flux_freq = pair[0][1:]
        flux_val = pair[1]
        is_flux = True
    else:
        is_flux = False
if is_flux == True:
    flux_density_str = 'Flux density at {} MHz: {} mJy\n'.format(flux_freq, flux_val)
else:
    flux_density_str = ''

pta_list = ''
for pta in pta_database.keys():
    if psr_name in pta_database[pta] or psr_bname in pta_database[pta]:
        pta_list += pta + ', '
    
output = 'Pulsar: {}\n'.format(psr_name) +\
         'RA: {}\n'.format(psr_ra) +\
         'Dec: {}\n'.format(psr_dec) +\
         'Period: {:.3e} s\n'.format(psr_period.value, 3) +\
         'Pdot: {:.3e}\n'.format(psr_pdot) +\
         'DM: {} pc/cm3\n'.format(psr_dm.value)
         
if psr_rm != None:
    output += 'RM: {} rad/m^2\n'.format(psr_rm)
    
output += flux_density_str +\
         'Characteristic age: {:.3e} yr\n'.format(psr_char_age.value) +\
         'Surface magnetic field: {:.3e} G'.format(round(psr_B_S.value, 3)) +\
         visible_telescopes

if pta_list != '':
    pta_list = '\nPTAs: ' + pta_list
    if pta_list[-1] == ' ':
        pta_list = pta_list[:-2]
    output += pta_list + '\n'

if len(output + wiki_link) > 280:
    print('Wikipedia link exists but will not be included; length of Tweet' +\
          ' would be {} characters.'.format(len(output + wiki_link)))
else:
    output += wiki_link
    
Output = 'TEST' + '\n' + output

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
#for i in range(len(pers)):
#    p = pers[i]
#    pdot = pdots[i]
#    if pdot > 10**(-16):
#        color = '#0057b7'
#    else:
#        color = '#ffd700'
#    plt.scatter([p], [pdot], c=color)
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
#plt.savefig(directory + '\\ppdot.png', bbox_inches="tight")
#plt.show()

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
#for i in range(len(l_list)):
    #l = l_list[i]
    #b = b_list[i]
    #if b > 0:
    #    color = '#0057b7'
    #else:
    #    color = '#ffd700'
    #ax2.scatter(l, b, c=color)
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
#plt.savefig(directory + '\\skymap.png', bbox_inches="tight")
#plt.show()

fig = plt.gcf()
fig.set_size_inches(12.0, 7.0)

def get_profile(jname, bname):
    is_webpage = False
    for name in [jname, bname]:
        url = 'http://rian.kharkov.ua/decameter/EPN/epndb/{}/'.format(name)
        response = requests.get(url)
        if response.status_code == 200:
            is_webpage = True
            
            # Courtesy of https://stackoverflow.com/a/24618186/6535830
            html = urlopen(url).read()
            soup = BeautifulSoup(html, features="html.parser")

            # kill all script and style elements
            for script in soup(["script", "style"]):
                script.extract()    # rip it out

            # get text
            text = soup.get_text()

            # break into lines and remove leading and trailing space on each
            lines = (line.strip() for line in text.splitlines())
            # break multi-headlines into a line each
            chunks = (phrase.strip() for line in lines for phrase in line.split("  "))
            # drop blank lines
            text = '\n'.join(chunk for chunk in chunks if chunk)
            
            text = text.split('\n')
            freqs = []
            for line in text:
                if 'GHz' in line:
                    freq = line.replace(' GHz', '')
                    freqs.append(float(freq))
            
            chosen_freq = random.choice(freqs)
            chosen_index = text.index('{} GHz'.format(chosen_freq))
            survey_index = chosen_index + 5
            survey_name = text[survey_index].split(' ')[1]
            
            mhz_freq = int(float(chosen_freq) * 1000)
            image_url = 'http://rian.kharkov.ua/decameter/EPN/epndb/{}/{}_{}.epn.gif'.format(name, survey_name, mhz_freq)
            with open('{}\\tweeted_pulsars\\{}_profile.gif'.format(directory, jname), 'wb') as f:
                f.write(requests.get(image_url).content)
            print('Image {}_profile.gif downloaded.'.format(jname))
            
    if is_webpage != True:
        print('No image available.')
        
    return is_webpage

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
        # Solution from https://stackoverflow.com/a/43660687/6535830
        filenames = ["{}\\tweeted_pulsars\\{}.png".format(directory, psr_name)]
        media_ids = []
        # The EPN website with profiles is based in Kharkiv and
        # is currently down due .
        #if get_profile(psr_name, psr_bname) == True:
        #    filenames.append("{}\\tweeted_pulsars\\{}_profile.gif".format(directory, psr_name))
        for filename in filenames:
            res = api.media_upload(filename)
            media_ids.append(res.media_id)
        api.update_status(status=output, media_ids=media_ids)
        print("Tweeting successful!")
    except tweepy.error.TweepError as e:
        print("Tweepy error!")
        print(e)
        #print("Invalid account credentials.")
