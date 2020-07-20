from __future__ import division, print_function
import argparse, sys

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

  source_altaz = source.transform_to(AltAz(location=obs_loc, obstime=obs_time))
  
  visible=False
  if source_altaz.alt.value > ob[3]:
    visible=True
  
  print("Source visible : {0}".format(visible))

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
  args = a.parse_args()
  main()

