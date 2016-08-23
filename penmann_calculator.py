#!/usr/bin/python
"""
Created on Sun Dec  6 10:12:46 2015

@author: theodor
"""

import math, datetime
import numpy as np

class pen_daily:
    '''
    calculate daily ET as it is described in the FAO paper Irrigation and
    drainage paper number 56 and displayed in box 11 example 18 
    http://www.kimberly.uidaho.edu/water/fao56/fao56.pdf
    '''
    def __init__(self,Tmax,Tmin,RHmax,RHmin,u_mean,Rs_mean,timestamp,alt,lat): 
        '''
        These are all daily parameters - max, min or mean as mentioned
        alt in [m]
        lat in decimal degrees
        '''
        self.Tmax   = float(Tmax)       # daily max temp (C)
        self.Tmin   = float(Tmin)       # daily min temp (C)
        self.RHmax  = float(RHmax)      # daily max RH(%)
        self.RHmin  = float(RHmin)      # daily min RH(%)
        self.u      = float(u_mean)     # daily mean wind speed
        self.Rs     = Rs_mean*0.0864         # short radiation [Wm-2] to [MJm-2*day-1] (daily mean)
        self.year   = timestamp.year    # year
        self.doy    = timestamp.timetuple().tm_yday # day of year
        self.month  = timestamp.month   # month
        self.alt    = float(alt)        # elevation of location(m)
        self.Gsc    = 0.0820            # solar constant [MJ m-2 min-1]
        self.sigma  = 4.903E-9  # Stefan-Bpltzman constant [MJK-4 m-2 day-1]
        self.phi    = math.radians(float(lat)) # decimal degrees to radians
        self.albedo = 0.23              # albedo
        self.G      = 0.0               # soil heat flux ~ 0 for daily ET 
    def t_mean(self): # [eq. 9]
        '''
        temp. mean of extremes
        '''
        return (self.Tmax+self.Tmin)/2.0
        
    def pressure(self): #[eq. 7]
        '''
        function that calculates pressure
        based on altitude returns pressure [kPa]
        '''
        p1 = (293-0.0065*self.alt)/293
        return 101.3*(p1**5.26)
 
    def gamma(self): # [eq.13]
        '''
        function that calculated gamma (psychrometric constant)
        based on pressure returns gamma
        '''
        return 0.665E-3*self.pressure() 
        
    def slope(self,T): # [eq. 13]
        '''
        the slope of the relationship
        between saturation vapour pressure and temperature
        '''
        p1   = 4098*(0.6108*math.exp(17.27*T/(T+237.3)))
        p2 = (T+237.3)**2
        return p1/p2
    
    def e_svp(self,T): # [eq. 11]
        ''' 
        function that calculates saturation vapour pressure
        based on temperature returns saturation vapour presure
        '''
        return 0.6108*math.exp(17.27*T/(T+237.3))
        
    def e_a(self): # [eq 17]
        '''
        Actual vapour pressure (ea) [kPa]
        derived from relative humidity data
        '''
        a1 = self.e_svp(self.Tmin)*self.RHmax/100
        a2 = self.e_svp(self.Tmax)*self.RHmin/100
        return (a1+a2)/2

    def e_s(self): # [eq. 12]
        '''
        Saturated vapour pressure (es) [kPa]
        '''
        a1 = self.e_svp(self.Tmin)
        a2 = self.e_svp(self.Tmax)
        return (a1+a2)/2
        
    def vp_defic(self):  
        '''
        vapor pressure deficit
        '''
        return self.e_s() - self.e_a()
        
    def dr(self): # [eq.23]
        ''''
        inverse relative distance Earth-Sun
        '''
        return 1 + 0.033*math.cos(2*math.pi*self.doy/365)
    
    def delta(self): # [eq. 24]
        '''
        solar declination [rad]
        '''      
        return 0.409*math.sin(2*math.pi*self.doy/365 - 1.39)    
        
    def omega(self): # [eq.25]
        '''
        sunset hour angle [rad] 
        '''
        delta = self.delta()
        phi = self.phi
        return math.acos(-math.tan(phi)*math.tan(delta))
       
    def rad_ext(self): # [eq. 21]
        '''
        extraterrestrial radiation (Ra) [MJ m-2 day-1]
        '''
        Gsc = self.Gsc
        dr  = self.dr()
        omega  =self.omega()
        phi = self.phi
        delta = self.delta()
        
        p1 = 24*60/math.pi
        p2 = Gsc* dr
        p3 = (omega*math.sin(phi)*math.sin(delta)) + (math.cos(phi)*math.cos(delta)*math.sin(omega))
        return p1*p2*p3
        
    def daylight_hours(self): #[eq. 34]
        '''
        daylight hours (N)
        '''
        omega = self.omega()
        return 24*omega/math.pi
        
    def rad_clear_sky(self): # [eq. 37]
        '''
        clear-sky radiation (Rso)
        '''
        Ra = self.rad_ext()
        return Ra*(0.75+2E-5*self.alt) 
        
    def rad_net_short(self): #[eq.38]
        '''
        Net Shortwave Radiation (Rns)
        '''
        return self.Rs*(1 - self.albedo) 
              
    def rad_net_long(self): # [eq. 39]
        Tmax = self.Tmax
        Tmin = self.Tmin
        Rs = self.Rs
        Rso = self.rad_clear_sky()
        e_a = self.e_a() 
        p1 = ((Tmax + 273.16)**4.0 + (Tmin + 273.16)**4.0)/2.0
        p2 = 0.34 - 0.14*math.sqrt(e_a)
        if Rs>Rso:
            p3 = 1
        else:
            p3 = 1.35*(Rs/Rso)-0.35
        return  self.sigma*p1*p2*p3
        
    def rad_net(self): # [eq.40]
        '''
        Net Radiation (Rn)
        '''
        Rns = self.rad_net_short()
        Rnl = self.rad_net_long()
        return Rns - Rnl
        
    def calc(self):
        '''
        function used to calculate the daily ETo
        returns ETo
        ET0 = [0.408 slope(Rn-G)+gamma*(900/(T+273.16))u2(vp_deficit)]/
              [slope + gamma(1+0.34u2)]
        '''  
        Rn = self.rad_net()
        T = self.t_mean()
        slope = self.slope(T)
        vp_def = self.vp_defic()
        u2 = self.u
        gamma = self.gamma()
        num = 0.408*slope*(Rn-self.G)+gamma*(900/(T+273.16))*u2*vp_def
        den = slope + gamma*(1+0.34*u2) 
        return num/den
        
#%%
#######################################################################################################################################        
class pen_hourly:
    '''
    calculate daily ET as it is described in the FAO paper Irrigation and
    drainage paper number 56 and displayed in example 19   
    http://www.kimberly.uidaho.edu/water/fao56/fao56.pdf
    '''
    def __init__(self, T, RH, u, Rs, timestamp, alt, lon, lat): 
        '''
        These are all hourly parameters - e.g: data for 14:00 is
        the mean or sum (for rain) of 14:00 - 14:59:59
        timestamp for winter time
        alt in [m]
        lat and lon in decimal degrees
        Lz and Lm written for Israel timezone
        '''
        self.T      = float(T)          # mean hourly temp (C)
        self.RH     = float(RH)         # mean hourly RH(%)
        self.u      = float(u)          # hourly mean wind speed at 2m
        self.Rs     = Rs*3.6E-3         # short radiation [Wm-2] to [MJm-2*h-1] (daily mean)
        self.year   = timestamp.year    # year
        self.doy    = timestamp.timetuple().tm_yday # day of year
        self.month  = timestamp.month   # month
        self.hour   = timestamp.hour    # hour
        self.alt    = float(alt)        # elevation of location(m)
        self.Gsc    = 0.0820            # solar constant [MJ m-2 min-1]
        self.sigma  = 2.043E-10         # Stefan-Bpltzman constant [MJK-4 m-2 h-1]
        self.phi    = math.radians(float(lat)) # decimal degrees to radians 
        self.albedo = 0.23              # albedo     
        self.Lz     = 330.0 # longitude of the centre of the local time zone (for Israel) [degrees west of Greenwich]
        self.Lm     = 360.0 - float(lon)  # longitude of the measurement site [degrees west of Greenwich], 
                                          # if Israel lon is decimal deg east of Greenwich 
#        self.Lm     = (lon)  # longitude of the measurement site [degrees west of Greenwich]

        return
        
    def dr(self): # [eq.23]
        ''''
        inverse relative distance Earth-Sun
        '''
        return 1 + 0.033*math.cos(2*math.pi*self.doy/365)
        
    def delta(self): # [eq. 24]
        '''
        solar declination [rad]
        '''      
        return 0.409*math.sin(2*math.pi*self.doy/365 - 1.39)  
    
    def omega_s(self): # [eq. 25]
        '''
        The sunset hour angle
        '''
        phi = self.phi 
        delta = self.delta()
        return math.acos(-math.tan(phi) * math.tan(delta))
        
    def daylight_hours(self): #[eq. 34]
        '''
        daylight hours (N)
        '''
        omegas = self.omega_s()
        return 24*omegas/math.pi   
        
    def omega(self): # [eq. 31]
        t = self.hour + 0.5
        Lz = self.Lz
        Lm = self.Lm
        Sc = self.Sc()
        return (math.pi/12)*((t + 0.06667*(Lz - Lm) + Sc) - 12)
        
    def Sc(self): # [eq. 32]
        '''
        seasonal correction for solar time
        '''
        b = (2*math.pi*(self.doy - 81))/364
        return 0.1645*math.sin(2*b) - 0.1255*math.cos(b) - 0.025*math.sin(b)

    def rad_clear_sky(self): # [eq. 37]
        '''
        clear-sky radiation (Rso)
        '''
        Ra = self.rad_ext()
        return Ra*(0.75+2E-5*self.alt) 
        
    def rad_ext(self): # [eq. 28]
        '''
        extraterrestrial radiation (Ra) [MJ m-2 h-1]
        '''
        Gsc = self.Gsc
        dr  = self.dr()
        omega1  = self.omega() - (math.pi*1/24.0) # 1 stands for 1 hour 
        omega2  = self.omega() + (math.pi*1/24.0) # 1 stands for 1 hour
        phi = self.phi
        delta = self.delta()
        
        p1 = 12*60/math.pi
        p2 = Gsc* dr
        p3 = (omega2 - omega1)*math.sin(phi)*math.sin(delta) + \
            math.cos(phi)*math.cos(delta)*(math.sin(omega2) - math.sin(omega1))
        return p1*p2*p3

    def rad_net_long(self): # [eq. 39 for hourly]
        '''
        Net Longwave Radiation (Rnl)
        '''
        T = self.T
        Rs = self.Rs
        Rso = self.rad_clear_sky()
        e_a = self.e_a() 
        p1 = (T + 273.16)**4.0 
        if e_a < 0.0:
            print 'e_a < 0'
        p2 = 0.34 - 0.14*math.sqrt(e_a)
        if Rs>Rso:
            p3 = 1
        else:
            p3 = 1.35*(Rs/Rso)-0.35
        return  self.sigma*p1*p2*p3       
        
    def rad_net_short(self): #[eq.38]
        '''
        Net Shortwave Radiation (Rns)
        '''
        return self.Rs*(1 - self.albedo) 

    def rad_net(self): # [eq.40]
        '''
        Net Radiation (Rn)
        '''
        Rns = self.rad_net_short()
        Rnl = self.rad_net_long()
        return Rns - Rnl

    def soil_heat_flux(self): # [eq. 45-46]
        '''
        soil heat flux
        '''
        Rn = self.rad_net()
        if Rn <= 0.0:
            G = Rn*0.5
        else:
            G = Rn*0.1          
        return G
        
    def e_svp(self,T): # [eq. 11]
        ''' 
        function that calculates saturation vapour pressure
        based on temperature returns saturation vapour presure
        '''
        return 0.6108*math.exp(17.27*T/(T+237.3))
        
    def e_a(self): # [eq 54]
        '''
        Actual vapour pressure (ea) [kPa]
        derived from relative humidity data
        '''
        return self.e_s()*self.RH/100
 
    def e_s(self): # [eq. 12]
        '''
        Saturated vapour pressure (es) [kPa]
        '''
        return self.e_svp(self.T)

    def pressure(self): #[eq. 7]
        '''
        function that calculates pressure
        based on altitude returns pressure [kPa]
        '''
        p1 = (293-0.0065*self.alt)/293
        return 101.3*(p1**5.26)
        
    def vp_defic(self):  
        '''
        vapor pressure deficit
        '''
        return self.e_s() - self.e_a() 


    def slope(self,T): # [eq. 13]
        '''
        the slope of the relationship
        between saturation vapour pressure and temperature
        '''
        p1   = 4098*(0.6108*math.exp(17.27*T/(T+237.3)))
        p2 = (T+237.3)**2
        return p1/p2

    def gamma(self): # [eq.13]
        '''
        function that calculated gamma (psychrometric constant)
        based on pressure returns gamma
        '''
        return 0.665E-3*self.pressure() 
        
    def calc(self): # [eq. 53]
        '''
        function used to calculate the daily ETo
        returns ETo
        ET0 = [0.408 slope(Rn-G)+gamma*(37/(T+273.16))u2(vp_deficit)]/
              [slope + gamma(1+0.34u2)]
        '''
        params = np.array([self.T,
                           self.RH,
                           self.u,
                           self.Rs])
        if  any(np.isnan(params)) == True:
            ET = np.nan
        else:
            Rn = self.rad_net()
            T = self.T
            slope = self.slope(T)
            vp_def = self.vp_defic()
            u2 = self.u
            G = self.soil_heat_flux()
            gamma = self.gamma()
            num = 0.408*slope*(Rn - G)+gamma*(37/(T + 273.16))*u2*vp_def
            den = slope + gamma*(1 + 0.34*u2)
            ET = num/den
            if ET < 0.0:
                ET =0.0
        return ET
        
#%%
def adjust_ws(u,z):
    '''
    adjusts wind speed to 2m height 
    u = unadjusted wind speed
    z = elevation
    '''
    u2 = u * 4.87/(math.log(67.8 * z - 5.42))
    return u2
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    