
"""
Created on Mon Jul 13 11:59:34 2020

@author: Sunstang Strategy Team
"""

import math
import os
import datetime
import random
import numpy as np
import pandas as pd 

def Aero_Power(A, Cd, p, V): #Aerodynamic power loss calculation
    return 0.5*p*(V/3.6)**3*A*Cd

def Roll_Resist(Crr,V,M): #Calculate the power loss from friction/rolling resistance
    return Crr*(1+V/161)*M*9.81*V/3.6
 
def Array_Power(Day, Latitude, Time, Sunrise, DayLength, Pmax, Driving): #Calculate the power gain from the solar array
    SLL = 23.5*math.sin(math.radians((180*(Day-82))/182.5))
    phi_N = Latitude - SLL
    phi = 90 - ((90-phi_N)*math.sin(math.radians(180*(Time-Sunrise)/DayLength)))
    if Driving == True:
        theta = phi
    
    return Pmax*(math.cos(math.radians(phi))**0.3)*math.cos(math.radians(theta))

def Grav_Power(V,M,d,alt1,alt2): #Calculate the power loss due to gravitatational effects
    return V*M*9.81*math.sin(math.atan((alt2-alt1)/1000/d))/3.6

def Haversine(lat1, lat2, long1, long2, r):
    lat1=math.radians(lat1)
    lat2=math.radians(lat2)
    long1=math.radians(long1)
    long2=math.radians(long2)
    return 2*r*math.asin(math.sqrt(pow(math.sin((lat2-lat1)/2),2)+(math.cos(lat1)*math.cos(lat2)*pow(math.sin((long2-long1)/2),2))))

def delta_T(speed,distance):        #speed = km/h, distance = km
    return datetime.timedelta(seconds = distance/speed*3600)

def Batt_Power(Pdrag,Prr,Pg,Pk,Parr,MotorEff,Pelec): 
    return ((Pdrag+Prr+Pg+Pk)/MotorEff)+Pelec-Parr

def Energy(Power, Time):
    return (Power*Time)/1000

def date_to_jd(year,month,day):
    # Convert a date to Julian Day.
    # Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet', 
    # 4th ed., Duffet-Smith and Zwart, 2011.
    # This function extracted from https://gist.github.com/jiffyclub/1294443
    if month == 1 or month == 2:
        yearp = year - 1
        monthp = month + 12
    else:
        yearp = year
        monthp = month
    # this checks where we are in relation to October 15, 1582, the beginning
    # of the Gregorian calendar.
    if ((year < 1582) or
        (year == 1582 and month < 10) or
        (year == 1582 and month == 10 and day < 15)):
        # before start of Gregorian calendar
        B = 0
    else:
        # after start of Gregorian calendar
        A = math.trunc(yearp / 100.)
        B = 2 - A + math.trunc(A / 4.)

    if yearp < 0:
        C = math.trunc((365.25 * yearp) - 0.75)
    else:
        C = math.trunc(365.25 * yearp)
    D = math.trunc(30.6001 * (monthp + 1))
    jd = B + C + D + day + 1720994.5
    return jd    

def sunrise(latitude, longitude, timezone, date):
    latitude = math.radians(latitude)
    longitude = math.radians(longitude)
    
    #constants
    jd2000 = 2451545 #the julian date for Jan 1 2000 at noon
    earth_tilt = math.radians(23.44)
    sun_disc = math.radians(-0.83)

    #equations
    jd_now = date_to_jd(date.year,date.month,date.day) #current julian date
    n = jd_now - jd2000 + 0.0008 #Current julian day
    jstar = n - math.degrees(longitude)/360 #Mean solar noon
    M = math.radians(math.fmod(357.5291 + 0.98560028 * jstar,360)) #Solar mean anomaly - degrees
    C = 1.9148 * math.sin(M) + 0.0200 * math.sin(2*M) + 0.0003 * math.sin(3*M) #Equation of the center
    lamda = math.radians(math.fmod(math.degrees(M) + C + 180 + 102.9372,360)) #Eliptic longitude - degrees
    Jtransit = 2451545.5 + jstar + 0.0053 * math.sin(M) - 0.0069 * math.sin(2*lamda) #Solar transit
    angle_delta = math.asin(math.sin(lamda) * math.sin(earth_tilt)) #Deinclination of the sun
    omega = math.acos((math.sin(sun_disc) - math.sin(latitude) * math.sin(angle_delta))/(math.cos(latitude) * math.cos(angle_delta))) #Hour angle

    Jrise = Jtransit - math.degrees(omega)/360
    Jset = Jtransit + math.degrees(omega)/360

    numdays_rise = Jrise - jd2000 +0.5 +timezone/24
    numdays_set = Jset - jd2000+0.5 +timezone/24
    sunrise = datetime.datetime(2000, 1, 1) + datetime.timedelta(numdays_rise)
    sunset = datetime.datetime(2000, 1, 1) + datetime.timedelta(numdays_set)
    global DL 
    DL = sunset-sunrise
    return sunrise

def conv_to_DT(H, M, S, MS):    #Convert time to military decimal time
    return H+(M/60)+(S/3600)+(MS/3600000000)

def Kine_Power(V_now, V_past, Dist, M, Time): 
    if Time.hour == 9 and Time.minute == 0 and Time.second == 0 and Time.microsecond == 0:
        V_past = 0
    return 5.46e-7*M*9.81*(((V_now**2-V_past**2)*(V_now+V_past))/(Dist))

def generateData(route_data_df,vehicle_specs_csv): #Returns a saved csv file
    
    Num_Segments = len(route_data_df-1) #calculate the number of race route segments
    
    #Object being returned
    emm_variables = pd.DataFrame(index = [i for i in range(Num_Segments)],columns = ["Segment Distance (km)","Segment Velocity (km/h)","Segment Start Time","Segment Elapsed Time (s)","Segment End Time", "Array Power (W)", "Rolling Power (W)", 
    "Gravitaional Power (W)",'Drag (W)',"Kinetic Power (W)", "Battery Power (W)", "Battery Energy (kWh)",'Energy Difference'])

    #Constants
    p = 1.17    #density of air
    Max_Array_Power = vehicle_specs_csv['Value']['Max_Array_Power']
    MotEff =  vehicle_specs_csv['Value']['MotEff']
    Loaded_weight = vehicle_specs_csv['Value']['Loaded_weight']
    A = vehicle_specs_csv['Value']['A']
    Cd = vehicle_specs_csv['Value']['Cd']
    Crr = vehicle_specs_csv['Value']['Crr']
    
    #Date & Time Variables
    Start_Day = '2021-10-22T09:00:00'                       #Start day & time for race in str, will eventually be a value read from a csv file
    Start_Time = datetime.datetime.fromisoformat(Start_Day)       #create datetime object from the Stat_Day str
    Timezone = 9.5    #Timezone of the race, will eventually be read in from a csv file
    Current_Time =  Start_Time  
    Batt_energy = 0    
    count = 0
   
    for x,row in route_data_df.iterrows():
        if x == len(emm_variables):
            break
        #Load in data
        SR = sunrise(route_data_df['latitude'][x], route_data_df['longitude'][x], Timezone, Current_Time.date())
        
        #Determine velocity, distance travelled and time taken
        if x == 0:
            emm_variables['Segment Start Time'][x] = Start_Time
        else: 
            emm_variables['Segment Start Time'][x] = emm_variables['Segment End Time'][x-1]
        
        emm_variables['Segment Velocity (km/h)'][x] = random.randrange(25,88)  #Solar Car Speed **This needs to be a list/array with a size of the number of data points along the route - 1]
        emm_variables['Segment Distance (km)'][x] = Haversine(route_data_df['latitude'][x],route_data_df['latitude'][x+1], route_data_df['longitude'][x], route_data_df['longitude'][x+1], 6371)
        dT = delta_T(emm_variables['Segment Velocity (km/h)'][x],emm_variables['Segment Distance (km)'][x])
        emm_variables['Segment Elapsed Time (s)'][x] = dT.total_seconds()
        emm_variables['Segment End Time'][x] = emm_variables['Segment Start Time'][x]+dT
    
        #Determine powers and energy
        emm_variables['Array Power (W)'][x] = Array_Power(Current_Time.timetuple().tm_yday,route_data_df['latitude'][x],conv_to_DT(Current_Time.hour,Current_Time.minute,Current_Time.second,Current_Time.microsecond),conv_to_DT(SR.hour,SR.minute,SR.second,SR.microsecond), DL.total_seconds()/3600,Max_Array_Power,1)
        emm_variables['Drag (W)'][x] = Aero_Power(A,Cd,p,emm_variables['Segment Velocity (km/h)'][x])
        emm_variables['Rolling Power (W)'][x] = Roll_Resist(Crr,emm_variables['Segment Velocity (km/h)'][x],Loaded_weight)
        emm_variables['Gravitaional Power (W)'][x] = Grav_Power(emm_variables['Segment Velocity (km/h)'][x], Loaded_weight, emm_variables['Segment Distance (km)'][x], route_data_df['altitude'][x], route_data_df['altitude'][x+1])
        
        if x == 0:
            emm_variables['Kinetic Power (W)'][x] = Kine_Power(emm_variables['Segment Velocity (km/h)'][x], 0, emm_variables['Segment Distance (km)'][x], Loaded_weight, Current_Time)
        else:
            emm_variables['Kinetic Power (W)'][x] = Kine_Power(emm_variables['Segment Velocity (km/h)'][x], emm_variables['Segment Velocity (km/h)'][x-1], emm_variables['Segment Distance (km)'][x], Loaded_weight, Current_Time)
        
        emm_variables['Battery Power (W)'][x] = Batt_Power( emm_variables['Drag (W)'][x],  emm_variables['Rolling Power (W)'][x],  emm_variables['Gravitaional Power (W)'][x],  emm_variables['Kinetic Power (W)'][x],  emm_variables['Array Power (W)'][x],MotEff, vehicle_specs_csv['Value']['Elec_consumption'])
        emm_variables['Battery Energy (kWh)'][x] = Energy(emm_variables['Battery Power (W)'][x], emm_variables['Segment Elapsed Time (s)'][x]/3600)
        
        Batt_energy += emm_variables['Battery Energy (kWh)'][x] 
        Batt_Energy_Ave = Batt_energy/(x+1)
        emm_variables['Energy Difference'][x] = abs(Batt_Energy_Ave-emm_variables['Battery Energy (kWh)'][x] )
        count+=1
        print(count)
        #Runs if a whole day has passed
        if emm_variables['Segment End Time'][x].hour == 18:
            Current_Time.replace(day=Current_Time.day+1, hour=9, minute=0, second=0, microsecond=0)
        else:
            Current_Time = emm_variables['Segment End Time'][x]
    return emm_variables

#Data directory where the files are stored
cwd = os.getcwd()
print("The current working directory is: " + cwd)
data_dir = cwd + '\sunstang_data'

#Load data sets in
route_data_filename = input("Please input which competition data you would like to access: ")
route_data_df = pd.read_csv(os.path.join(data_dir,route_data_filename))

#Vehicle Specifications
vehicle_specs_filename = input("Please input vehicle specification file name: ")
vehicle_specs_csv = pd.read_csv(os.path.join(data_dir,vehicle_specs_filename),index_col=0)
print("Specs and data sets have been uploaded")                             

Req_Datasets = int(input("Required number of datasets: "))
count = 0
data = pd.DataFrame()
while count != Req_Datasets:
    data = generateData(route_data_df, vehicle_specs_csv)
    count +=1
    data.to_csv(os.path.join(data_dir,"Generated Dataset Number " + str(count) + '.csv'))
    print('Generated dataset number ' + str(count))
    