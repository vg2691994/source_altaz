from __future__ import division, print_function
import argparse, sys
import numpy as np
from astropy.coordinates import AltAz, Angle, EarthLocation, SkyCoord, FK5
from astropy.time import Time
from astropy import units as u

from parse_cfg import parse_cfg as pcfg

def get_source(name, source_list):
  if name in source_list.keys():
    ra = Angle(source_list[name][0], unit=u.hourangle)
    dec = Angle(source_list[name][1], unit=u.degree)
    
    source = SkyCoord(ra, dec, frame=FK5(equinox="J2000"))
    return source
  else:
    return None

def main():
  obs_locations = pcfg("observatories.cfg", delim="\s+", return_dict=True)
  ob = obs_locations[args.obs_loc]
  obs_loc = EarthLocation.from_geocentric(ob[0],ob[1], ob[2], unit=u.meter)

  source_list = pcfg("sources.list", delim="\s+", return_dict=True)
  if args.source_name is not None:
    source = get_source(args.source_name, source_list)
    if source is None:
      try:
        source = SkyCoord.from_name(args.source_name)
      except:
        raise ValueError("Could not find source name {0} in {1} or in astropy's SIMBAD databse".format(args.source_name, "sources.list"))

  elif args.source_coords is not None:
    ra = Angle(args.source_coords[0], unit=u.hourangle)
    dec = Angle(args.source_coords[1], unit=u.degree)
    source = SkyCoord(ra, dec, frame=FK5(equinox="J2000"))

  obs_time = Time(args.obs_time) 

  #Get the alt_az values for the whole LST day, centered around obs_time
  seconds_per_sidereal_day = 86164.09
  times = obs_time + (np.arange(0, seconds_per_sidereal_day, args.precision) - seconds_per_sidereal_day/2) * u.second
  
  source_altaz = source.transform_to(AltAz(location=obs_loc, obstime=obs_time))
  all_alt_az = source.transform_to(AltAz(location=obs_loc, obstime=times)) 
  
  #Find the times when the source rises and sets

  signs = np.sign(all_alt_az.alt.value - ob[3])
  diffs = np.diff(signs)

  max_alt_idx = np.argmax(all_alt_az.alt.value)
  max_alt = "{0:.1f} degrees".format(all_alt_az.alt.value[max_alt_idx])
  max_alt_time = str(times[max_alt_idx])

  if not np.any(diffs):
    #All values inside signs must be same, therefore, the source never crossed the horizon

    if np.all(signs == 1):
      visible = True
      rise_time = "Always up"
      set_time = "Always up"
    
    elif np.all(signs == -1):
      visible = False
      rise_time = "Never up"
      set_time = "Never up"

  else:
    #Source would have crossed the horizon twice. The diffs array would contain all 0s except for two places which should have +2 and -2.
    rise_idx = np.where(np.diff(np.sign(all_alt_az.alt.value - ob[3])) > 0)[0][0]
    set_idx = np.where(np.diff(np.sign(all_alt_az.alt.value - ob[3])) < 0)[0][0]

    rise_time = str(times[rise_idx])
    set_time = str(times[set_idx])

    visible=False
    if source_altaz.alt.value > ob[3]:
      visible=True

  print("At the requested obs_time(UTC) - {0}:".format(obs_time))
  print("Source visible : {0}".format(visible))
  print("Source altitude : {0}".format(source_altaz.alt))
  print("Source azimuth : {0}".format(source_altaz.az))

  print("\n")
  print("Nearest rise time(UTC): {0}".format(rise_time))
  print("Nearest set time(UTC): {0}".format(set_time))
  print("Max_altitide acheived: {0} at UTC: {1}".format(max_alt, max_alt_time))

if __name__ == '__main__':
  for i, arg in enumerate(sys.argv):
    if (arg[0] == "-") and arg[1].isdigit():
      sys.argv[i] = ' ' + arg
  a = argparse.ArgumentParser()
  a.add_argument("-obs_loc", type=str, help="Name of the observatory", required=True)
  a.add_argument("-source_name", type=str, help="Name of source --\
      Either present in the list 'sources.list' or\
      Resolvable by astropy (uses SIMBAD)")
  a.add_argument("-source_coords", type=str, nargs=2, help="RA (hourangle) and Dec (degree) of the source")
  a.add_argument("-obs_time", type=str, help="Time of observation in UTC (fmt=YYYY-mm-ddTHH:MM:SS). Def = Time.now()", default = str(Time.now()))
  a.add_argument("-precision", type=float, help="Precision (in seconds) when computing the rise and set times (def = 60 secs)", default=60)
  args = a.parse_args()
  main()

