# script_name.py
#
# Usage: python3 script_name.py arg1 arg2 ...
#  Text explaining script usage
# Parameters:
#  arg1: description of argument 1
#  arg2: description of argument 2
#  ...
# Output:
#  A description of the script output
#
# Written by Brad Denby
# Other contributors: Allison Hai
#
# Optional license statement, e.g., See the LICENSE file for the license.

# import Python modules
import math # math module
import sys # argv

# "constants"
R_E_KM = 6378.137
E_E = 0.081819221456
# helper functions

## function description
def calc_denom(ecc, lat_rad):
  return math.sqrt(1.0-(ecc**2)*(math.sin(lat_rad)**2))

# initialize script arguments
o_x_km    = float('nan') #
o_y_km    = float('nan') #
o_z_km    = float('nan') #
x_km      = float('nan') #x coordinate 
y_km      = float('nan') #y coordinate  
z_km      = float('nan') #z coordinate 

# parse script arguments
if len(sys.argv)==7:
    o_x_km = float(sys.argv[1])
    o_y_km = float(sys.argv[2])
    o_z_km = float(sys.argv[3])
    x_km   = float(sys.argv[4])
    y_km   = float(sys.argv[5])
    z_km   = float(sys.argv[6])
else:
    print(\
        'Usage: '\
            'python3 o_x_km o_y_km o_z_km z_km y_km z_km ...'\
    )
    exit()

# determine the ECEF vector from the station to the object
#ECEF vector from the station to the object
ecef_v  = [x_km-o_x_km, y_km-o_y_km,z_km-o_z_km] 
v1 , v2 , v3 = ecef_v 


## determine the ground site latitude, longitude, and height above the ellipsoid 
# calculate longitude
lon_rad = math.atan2(o_y_km,o_x_km)
lon_deg = lon_rad*180.0/math.pi

# initialize lat_rad, r_lon_km, r_z_km
lat_rad = math.asin(o_z_km/math.sqrt(o_x_km**2+o_y_km**2+o_z_km**2))
r_lon_km = math.sqrt(o_x_km**2+o_y_km**2)
prev_lat_rad = float('nan')

# iteratively find latitude
c_E = float('nan')
count = 0
while (math.isnan(prev_lat_rad) or abs(lat_rad-prev_lat_rad)>10e-7) and count<5:
  denom = calc_denom(E_E,lat_rad)
  c_E = R_E_KM/denom
  prev_lat_rad = lat_rad
  lat_rad = math.atan((o_z_km+c_E*(E_E**2)*math.sin(lat_rad))/r_lon_km)
  count = count+1
  
# calculate hae
hae_km = r_lon_km/math.cos(lat_rad)-c_E

#apply inverse rotations 
s_lat = math.sin(lat_rad) #cosine of latitude 
c_lat = math.cos(lat_rad)
s_lon = math.sin(lon_rad) #sin of longitude 
c_lon = math.cos(lon_rad)
#part 1, inverse R_z * inverse R_y
r_sez1 = [[(s_lat*c_lon) , (s_lat*s_lon + 0 +0) , (0 + 0 -c_lat)],    
         [(0  - s_lon +0) , (0 + c_lon + 0) , (0 + 0 + 0)],
         [(c_lat*c_lon + 0 + 0) , (c_lat*s_lon + 0 + 0) , (0 + 0 + s_lat)]] 
[[x1,x2,x3],[y1,y2,y3],[z1,z2,z3]]= r_sez1 
#part 2, multiply r_sez1 by ecef vector 
r_sez2 = [v1*x1+v2*x2+v3*x3, v1*y1+v2*y2+v3*y3, v1*z1+v2*z2+v3*z3]
s_km, e_km, z_km = r_sez2 
print(s_km)
print(e_km)
print(z_km)

