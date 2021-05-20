# PulsarOfTheDay

This repository houses the code for [the Twitter bot of the same name](https://twitter.com/PulsarOfTheDay). It's designed to randomly pick a pulsar from [the ATNF pulsar catalog](https://www.atnf.csiro.au/research/pulsar/psrcat/) and tweet out information about it - period, dispersion measure, characteristic age, etc. The bot also generates a couple of plots: a P-Pdot diagram and a skymap in galactic coordinates, both containing the chosen pulsar and a random selection of other pulsars.

The script is run on my computer via a .bat file, which did to one or two design choices, but you run it as you want from the command line. My copy of this script includes the Twitter account credentials, which I've obviously removed from this version. The default behavior is to Tweet, but you can just generate and display the output by using the `-local` flag, i.e.
```
$ python main.py -local
```
Without this flag, the script assumes you're trying to Tweet, and it will tell you it doesn't have the credentials, so make sure you use the `-local` flag if you're not me! There's also a second flag, `-manual`, which you can use to specify which pulsar the script should use. Simply type the pulsar's name (as it appears in the catalog) after the flag. For example:
```
$ python main.py -local -manual J0002+6216
```
Before running this, make sure you add in two lines: one specifying the directory the script and file are in, and one specifying the name of the catalog file. The database copy I have here is named `psrcat.db`, but if you download directly from ATNF's website, it *might* be named differently. In my script, then, I have the lines
```
directory = \path\to\file
database = 'psrcat.db'
```
The two non-standard libraries used here are Astropy and Tweepy; the latter really only matters if you're tweeting.
