# -*- coding: utf-8 -*-
"""
Download data from NOAA's "climate data online" repository, using API v2.
Information on this API is available here:
    https://www.ncdc.noaa.gov/cdo-web/webservices/v2

In order to use this script, you will have to request an access token from 
NOAA, and add it to the code at line 46.  Instructions are availalbe at the 
API link just provided.

Documentation on specific datasets is available at
https://www.ncdc.noaa.gov/cdo-web/datasets

Hourly data is not available through this API.  Information on hourly data can
be accessed here (International Surface Data): 
https://data.nodc.noaa.gov/cgi-bin/iso?id=gov.noaa.ncdc:C00532
https://www.ncdc.noaa.gov/isd

https://www.ncdc.noaa.gov/data-access/land-based-station-data/land-based-datasets

It looks like the hourly data is not available by API at all, but only by 
various manual searches.

Info on the Portland Station is here:
https://www.ncdc.noaa.gov/cdo-web/datasets/LCD/stations/WBAN:14764/detail (Portland)
https://www.ncei.noaa.gov/access/search/dataset-search?what=observationTypes:Land%20Surface%7CLand%2520Surface&placeId=3078&type=Zip%20Code&bbox=43.726,-70.336,43.661,-70.245&name=04103%252C%2520Portland%252C%2520ME%252C%2520USA&sdate=2018-01-01

Integrated Surface Hourly (ISH) -- apprently exposed as ISD

"""

import csv
import os
from os import path
from datetime import date, datetime
from time import sleep

import tkinter as tk
import tkinter.filedialog as filedialog
import tkinter.messagebox as msgbox

import requests


#Constants and Lookup lists
MYTOKEN = 'REPLACE THIS VALUE WITH AN ACCESS TOKEN'
MONTHLENGTHS = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
PORTLAND = 'USW00014764'
DAILY = 'GHCND:'
BASEURL = r'https://www.ncdc.noaa.gov/cdo-web/api/v2/'
DEFAULTFN = 'outfile.csv'

#TODO: Consider if the following should be selectable
# ENDPOINTS = [
#             'datasets' , #A  dataset is the primary grouping for data at NCDC.
#             'datacategories' , # A data category is a general type of data used to group similar data types.
#             'datatypes' , # A data type is a specific type of data that is often unique to a dataset.
#             'locationcategories' , # A location category is a grouping of similar locations.
#             'locations' , # A location is a geopolitical entity.
#             'stations' , # A station is a any weather observing platform where data is recorded.
#             'data' , # A datum is an observed value along with any ancillary attributes at a specific place and time.
#         ]

DATASETS = [
           'GHCND',   # Daily summaries
           'GSOM',   #Monthly Summaries
           'GSOY' #,   # Annual summaries
          # 'NORMAL_ANN',
          # 'NORMAL_DLY',
          # 'NORMAL_HLY',
          # 'NORMAL_MLY',
          # 'PRECIP_HLY',
          # 'PRECIP_15'
        ]


def getdatacategories():
    ''' Fetch listing of all available data categories.  Currently uncalled
    through the GUI, but included here as a useful utility.'''
    url = BASEURL + 'datacategories'
    parms = {'limit':'50'}
    headers = {'token': MYTOKEN}
    response = requests.get(url, headers = headers, params = parms)
    if not (response.status_code > 199 and response.status_code<300):
        print ('Server Error. Status Code:', response.status_code)
    else:
        return response.json()
    
def getdatasets():
    ''' Fetch listing of all available data sets.  Currently uncalled
    through the GUI, but included here as a useful utility.'''
    url = BASEURL + 'datasets'
    parms = {'limit':'50'}
    headers = {'token': MYTOKEN}
    response = requests.get(url, headers = headers, params = parms)
    if not (response.status_code > 199 and response.status_code<300):
        print ('Server Error. Status Code:', response.status_code)
    else:
        return response.json()  

def checkstation(code='GHCND:USW00014764'):
    ''' Check on the status of a specific station code'''
    #TODO:  Separate station code and data type code
    #TODO:  Add code to allow user to check stations for other data types
    # Need to figure out what those codes look like...
    url = BASEURL + 'stations' + r'/' + code
    headers = {'token': MYTOKEN}
    response = requests.get(url, headers = headers)
    if not (response.status_code > 199 and response.status_code<300):
        print ('Server Error. Status Code:', response.status_code)
    else:
        return response.json()   


def MaineStations():
    ''' Generate a list of all stations in Maine.'''
    url = BASEURL + 'stations' + '?locationid=FIPS:23&limit=500'
    headers = {'token': MYTOKEN}
    response = requests.get(url, headers = headers)

    if not (response.status_code > 199 and response.status_code<300):
        print ('Server Error. Status Code:', response.status_code)
    else:
        r = response.json()['results']
        return [(st['name'], st['id']) for st in r ]

def getmonthrange (year,month):
    monthlengths = MONTHLENGTHS
    start = date(year, month, 1).strftime('%Y-%m-%d')
    if (year/100.0 != year/100) and (year/4.0 == year/4):
        # it's a leapyear
        monthlengths[1]=29
    else:
        monthlengths[1]=28
    end = date(year, month, monthlengths[month-1]).strftime('%Y-%m-%d')
    return start, end

#Structure of the response:
#   Top level dictionary, consisting of "results" and "metadata" components
#   Metadata is another dictionary, containing a single entry, "resultset",
#       also a dict.
#   Results is a list of dictionaries, each containing:
            # 'datatype'
            # 'date'
            # 'station'
            # 'value'
            # 'attributes'
#   Attributes here appear to be codes for data quality, and they are partially
# described in the documentation.
            
# A data dictionary for the parameter codes I see at Portland
#        AWND = Average wind speed (m/s or mph) 
#        PRCP = Precipitation (mm or inches as per user preference) 
#        SNOW = Snowfall (mm or inches as per user preference) 
#        SNWD = Snow depth (mm or inches as per user preference) 
#        TAVG = Average  temperature  (Fahrenheit  or  Celsius)
#        TMAX = Maximum  temperature  (Fahrenheit  or  Celsius)
#        TMIN = Minimum  temperature
#        WDF2 = Wind direction fastest 2 minutes
#        WDF5 = Wind direction fastest 5 seconds
#        WSF2 = Wind speed fastest 2 minutes
#        WSF5 = Wind speed fastest 5 seconds

#%%
class MsgPanel(tk.Frame):
    ''' A utility class for a small message widget as part of the UI'''
    def __init__(self, master, row, column):
        #  constructor function builds the 
        tk.Frame.__init__(self, master,)   # Constructor for the parent class.
        self.config(relief = tk.GROOVE, bd = 1)
        tk.Label(self, text = 'Results:').grid(column=0, row=0,
                 sticky=(tk.W)) 
        self.msg = tk.Entry(self, relief=tk.SUNKEN,
                                  bd=5, width = 80)
        self.msg.grid(column=1, row=0, sticky=(tk.E, tk.W),
                            padx = 35, pady = 15)
        self.msg.insert(0,'Waiting for Parameters')
        
        self.grid(column=column, row=row, sticky=(tk.N, tk.W, tk.E, tk.S))

    def write(self, message):
        self.msg.delete(0, tk.END)
        self.msg.insert(0, message)
        self.update_idletasks()
        
    def clear(self, message):
        self.msg.delete(0, tk.END)
        self.update_idletasks()

#%%
class NOAAWeatherApp(tk.Frame):
    ''' Application to retreive weather data from NOAA's Climate Data Online
    (CDO) using API V2.  This provides access to daily weather summaries and
    climate normals, but not access to hourly weather data.
    
    Inherits from a Tkinter FRAME.
    
    Constructor requires one parameter, a  Tkinter instance that holds the
    NOAAWeatherAPP.
    
    This simplified application does not allow other (named) option parameters
    to be passed to the constructor, even though that is possible in tkinter.
    
    Class Methods:
        __init__             Initializes the App
        _buildtoolbar        Builds the toolbar which contains buttons
        _builddatapanel      Builds the panel that gathers the parameters for
                             the request.
        _buildmessagepanel   The message panel is the part of the app that
                             allows feedback to the user.
        _formatdate          This method is a callback to condition date entry
                             so that the resulting string is a well-formed date
        go                   Callback that fires off the request for data
        orderdata            Organizes one or many requests, month by month;
                             returns list of dicts containingdatas in long form
        retreivedata         Actually sends a request, returns dict derived
                             from JSON.  Data is in return['results']
                             as a list of dicts. Metadata in return['metadata']

        appendtalldata       saves tall form data one data block at a time
        appendwidedata       saves wide form data one data block at a time.
        
    Class Attributes
        station             StringVar containing station code
        datasetlstbx        Listbox for dataset selection
        datform             StringVar indicating which data format you want
        sdate               StringVar for start date
        edate               StringVar for end date
        filename            StringVar for filename
        themsg              MsgPanel object for comunicating with the user
        firstdata           Flag to indicate if we are starting a new data file
    '''
    # Application inherits from Frame.  This version requires a toplevel master
    # here "root = Tk".
    # Most of the work is in creation of callback functions.
    # These are implemented as methods and properties of the App object.

    def __init__(self, master):
        # App constructor function builds the user interface

        tk.Frame.__init__(self, master)   # Constructor for the parent class.
        self.master.title('Online Access for NOAA Historic Weather Data ONLY WORKS FOR HOURLY')
        self.grid(column=0, row=0, sticky=(tk.N, tk.W, tk.E, tk.S))
        self.master.resizable(False, False)
        
        self._buildtoolbar()
        self._builddatapanel()
        self.themsg = MsgPanel(self, 5,0)

        
        # Create some padding
        for child in self.winfo_children():
            child.grid_configure(padx=5, pady=5)

    #######################
    # Functions to assemble the UI
    #######################
    def _buildtoolbar(self, row=0, column=0):
        '''Build a toolbar containing buttons'''
        toolbar = tk.Frame(self, relief = tk.GROOVE, bd = 1)
        # Build a couple of buttons in the toolbar
        quitbutton = tk.Button(toolbar, text="Quit", command=self.Quit)
        quitbutton.grid(column=0, row=0, padx=2, pady=2,
                               sticky=(tk.W))
        gobutton = tk.Button(toolbar, text="Go", command=self.go)
        gobutton.grid(column=1, row=0, padx=2, pady=2,
                               sticky=(tk.W))
        toolbar.grid(column=column, row=row, sticky=(tk.N, tk.W, tk.E, tk.S))
        
    def _builddatapanel(self, row=1, column=0):
        '''Build a Panel for Data Input, principally Station ID, dates, and
        filenames'''

        panel = tk.Frame(self, relief = tk.GROOVE, bd = 1)
        
        # Data request information first    
        # Station first  
        tk.Label(panel, text = 'Station', justify = 'left').grid(row = 0,
                column = 0, sticky = tk.NW)
        
        self.station = tk.StringVar() 
        self.station.set(DAILY+PORTLAND)
        stationentry = tk.Entry(panel, width = 25, textvar=self.station)
        stationentry.grid(row = 0, column = 1, sticky = tk.NW)
 
    
        # Now, which dataset do we want?  This is not connected to request yet...
        tk.Label(panel, text = 'Dataset', justify = 'left').grid(row = 1,
                column = 0, sticky = tk.NW)
        #TODO:  Convert this to a stringvar, for consistency (HOW?)
        self.datasetlstbx = tk.Listbox(panel, selectmode=tk.SINGLE,
                                 exportselection=False, height=6)
        self.datasetlstbx.grid(row = 1, column = 1)
        for item in DATASETS:
            self.datasetlstbx.insert(tk.END, item)
        self.datasetlstbx.select_set(0) # sets the first element, here 'GHCND'
    

        
        # Select data format to save  with radiobuttons
        self.datform = tk.StringVar()   # 
        self.datform.set('Wide Form')
        btnframe = tk.Frame(panel, relief = tk.GROOVE, bd=2)
        btnframe.grid(row = 0, column = 2, rowspan=2, sticky = tk.NE)
        tk.Label(btnframe, text ='Data Format').grid(row=0, column=0)
        tk.Label(btnframe, text = 'Long Form').grid(row=1, column=0)
        b1 = tk.Radiobutton(btnframe, variable = self.datform,
                            value = 'Long Form')
        b1.grid(row=1, column=1)
        tk.Label(btnframe, text = 'Wide Form').grid(row=2, column=0)
        b2 = tk.Radiobutton(btnframe, variable = self.datform,
                            value = 'Wide Form')
        b2.grid(row=2, column=1)

        #Calculate default date -- today
        now = datetime.now() # current date and time
        today = now.strftime('%Y-%m-%d')
        
        # Assemble dates for the request
        tk.Label(panel, text = 'Starting Date (YYYY-MM-DD):').grid(row = 2,
                column = 0, sticky = tk.W)
        self.sdate = tk.StringVar()
        self.sdate.set(today)
        sdateentry = tk.Entry(panel, textvariable=self.sdate, width = 25)
        sdateentry.bind('<KeyPress>', self._format_date)
        sdateentry.grid(row = 2, column = 1, sticky = tk.W)
        
        tk.Label(panel, text = 'Ending Date (YYYY-MM-DD):').grid(row = 2,
                column = 2, sticky = tk.W)
        self.edate = tk.StringVar()
        self.edate.set(today)
        edateentry = tk.Entry(panel, textvariable=self.edate, width = 25)
        edateentry.bind('<KeyPress>', self._format_date)
        edateentry.grid(row = 2, column = 3, sticky = tk.W)


        # Accept a filename
        
 
        tk.Label(panel, text = 'Output File Name:').grid(row = 4,
                column = 0, sticky = tk.W)
        self.filename = tk.StringVar()
        self.filename.set(DEFAULTFN)
        filenameentry = tk.Entry(panel, textvariable=self.filename, width = 82)
        filenameentry.grid(row = 4, column = 1, columnspan = 3, sticky = tk.W)
        fn_btn = tk.Button(panel, text='File Dialog', relief=tk.RAISED,
                        command=self._getfilename)
        fn_btn.grid(row = 4, column = 5)

        panel.grid(column=0, row=1, sticky=(tk.N, tk.W, tk.E, tk.S))   
        for child in panel.winfo_children():
            child.grid_configure(padx=5, pady=5)
    
    #######################
    # Callbacks of one sort or another.
    #######################
    
    def _format_date(self, event):
        '''Callback that restricts keystrokes as it's being entered.
        Allows entries that mach the following format:  'YYYY-MM-DD'
        '''
        
       # Works MOST of the time.  However, if hte user selects multiple chars
       # with the mouse, can get off-kilter.  The work around is to retype.....
       # Apparently, the alternate way to handle that would be to capture the
       # MOuse up event and remove the selection...

        if event.keysym in {'Tab', 'Return'} :   # allow to escape field
            return     #Allow
        
        widget = event.widget       # get the entry widget
        entry = widget.get()        # get the text content.  Note not via StringVar
        idx = widget.index(tk.INSERT)  # current cursor position
                
        if event.keysym == 'Left':  # Handle left arrows to skip dashes
            if idx == 5:
                widget.icursor(4)
            elif idx == 8:
                widget.icursor(7)
        
        elif event.keysym == 'Right':   # Handle right arrows to skip dashes
            if idx == 3:
                widget.icursor(4)
            elif idx == 6:
                widget.icursor(7)
        
        elif event.keysym == 'BackSpace':
            # convert a backspace to a 'Left' event, setting up to overwrite
            # at the correct spot
            widget.event_generate('<Left>')
            return 'break'
        
        elif event.keysym == 'Delete':
            # Silently ignore a delete
            # deletes tend to mess up  the order of characters
            # because this fxn overwrites values, maintaining positions
            return 'break'

        else:
            # block if char not a digit 
            # or date would be to long
            if not event.char.isdigit() or idx > 9:
                widget.bell()
                return 'break'  # block edit
            
            # Constrain to 19th, 20th, and 21st centuries
            # First digit must be 1 or 2
            if idx == 0:
                if event.char not in '12':
                    widget.bell()
                    return 'break'  # block edit
            if idx == 1:
                if entry[0] == '1' and not event.char in'89':
                    widget.bell()
                    return 'break'  # block edit
                if entry[0] =='2' and event.char !='0':
                    widget.bell()
                    return 'break'  # block edit

            # Constrain to valid months
            if idx == 5 or idx == 4:
                 if not event.char in '01':
                    widget.bell()
                    return 'break'  # block 
            if idx == 6:
                if entry[5] == '0' and event.char=='0':
                    widget.bell()
                    return 'break'  # block 
                if entry[5]=='1' and event.char not in '012':
                    widget.bell()
                    return 'break'  # block 
                
            # Constrain to Valid days
            if idx == 8 or idx == 7:
                if idx == 7:
                    widget.icursor(8)
                if event.char not in '0123':
                    widget.bell()
                    return 'break'  # block
            if idx == 9:
                if entry[8] == '0' and event.char=='0':
                    widget.bell()
                    return 'break'  # block 
                if entry[8] == '3' and event.char not in '01':
                    widget.bell()
                    return 'break'  # block
                
            # Finally, We get to ADDING or REPLACING new digits
            
            # ADDING. This shouldn't get called, since we pre-populated the
            # date field with the current date. It's here for completeness.
            if idx == widget.index(tk.END):
                # allow adding a new digit, inserting dashes where necessary             
                if idx == 4 or idx == 7:   # insert format char
                    widget.insert(idx,'-')
                    widget.icursor(idx+1)      # Move cursor before edit
    
            else:  # REPLACING a digit
                if idx in [4,7]:
                    # if overwriting dashes, skip to next character
                    widget.delete(idx, idx+2)  # Delete two
                    widget.insert(idx, '-')    # Add dash
                    widget.icursor(idx+1)      # Move cursor before edit
                else:  # ok to replace, so delete next char
                        widget.delete(idx)


    def _getfilename(self):
        '''Callback for selecting a file name and directory to save the data'''
        # Begin by defining options and opening the file dialog
        options = {}
        
        existingpath = self.filename.get()
        # If a file name other than the default already exists, use that as the
        # starting point
        if existingpath != DEFAULTFN:
            exdir, exfn = path.split(existingpath)
        else:
            exdir, exfn = os.getcwd(), DEFAULTFN
            # In python2, getcwd returns a string with double
            # backslashes. We need to convert for consistency with file dialog.
            exdir = exdir.replace('\\','/')

        options['initialdir'] = exdir
        options['initialfile'] = exfn
        options['title'] = 'Where to Save Downloaded Data?'
        
        #Run the dialog
        newpath = filedialog.asksaveasfilename(**options)
        if newpath:
            # filedialog returns a unicode string with single slashes in
            # the path. Need to convert to string.  This may fail if unicode that
            # is illegal in str is used. 
            newpath = str(newpath)
    
            # For visual simplicity, if we're in teh current dirtectory, just use
            # the filename
            newdir, newfile = path.split(newpath)
            if newdir == os.getcwd().replace('\\','/'):
                newpath = newfile
            
            #Set the filename
            self.filename.set(newpath)

           
# This function is just a spot to process requested data periods and send off
# requests.  Note that the structure here assumes we are requesting daily data.
# if not, there is no reason to call data mmonth by month.  For monthly and
# annual data (normals, suummaries), that's overkill.
    def go(self):
        '''
        Send off requests for data, and save, as specified in the GUI.
        
        The method calls self.orderdata, which pulls parameters from the app
        object (often the UI), ooks at the time period requested, and submits
        sequential HTPP requests to NOAA and assembles a single list of dicts.
        We read data by month to stay below normal 1000 record limit  and to
        reduce load on NOAA servers.
        '''
        self.themsg.write('Submitting Requests....')
        self.firstdata = True                          # A flag for first data saved       
        self.orderdata()       
        

    ######################################
    # Methods to retreive data from the NOAA Server and save them as CSVs
    ######################################
    def retreivedata(self,  start, end,
                     limit = 500, code = 'GHCND',   # Global Daily summaries
                     #units = 'metric'  # Caused problems....
                     ):
        '''Send single data request to NOAA and return data as a DICT derived
        from the JSON object returned by the server.
        
        This Dict includes both 'results' and 'metadata'. 'results' is a list
        of dicts contianing "long Form" data.
        '''     
        self.themsg.write('Request from %s to %s' % (start,end))
        payload = {'datasetid' : code,
                   'stationid': self.station.get(),
                   'startdate': start,
                   'enddate'  : end,
                   'limit'    : str(limit),   # max is 1000 -- roughly 60 days data
                   #'units'    : 'metric'  The units parameter causes errors?
                   }
        headers = {'token': MYTOKEN}
        response = requests.get(BASEURL+'data', headers = headers,
                                params = payload)
        #print(response.request.url)
        tries = 0
        finished = False
        problem=False    #TODO:  refactor to use python's exceptions
        while tries < 3 and not finished:
            tries += 1
            response = requests.get(BASEURL+'data', headers = headers,
                                params = payload)
            print('Requesting data begining %s' % start)
            if response.status_code == 200:
                finished = True
            else:
                problem=True
                print(response.status_code)
                self.themsg.write('Status Code: %s' % response.status_code)
                finished = True
                if response.status_code == 503:  # Service unavailable -- retry
                    self.themsg.write('Retrying.... {}'.format(tries+1))
                    finished = False
                    problem=False
                if response.status_code == 429:  # too many requests -- slow down
                    self.themsg.write('Retrying.... {}'.format(tries+1))
                    sleep(1)
                    finished = False
                    problem=False
            if tries==3:
                problem = True
            self.update_idletasks()
        if not problem and response.text:
            data= response.json()  # Extract response body and converts to dict
            return data
        else:
            return None
#%%
    def orderdata(self):
        ''' Dispatches orders for data in suitable intervals (usually monthly)
        assembles a data structure containing data from the start date to the
        end dat rettreived from the UI,. Optionally, and by default) saves
        the data to disk after each HTTP request is successful.
        '''
        
        # First, set up month range needed
        sd =datetime.strptime(self.sdate.get(), '%Y-%m-%d')
        ed =datetime.strptime(self.edate.get(), '%Y-%m-%d')
        
        if sd>ed:
            sd,ed = ed,sd
            self.themsg.write('Dates were reversed....  Corrected.')
            
        sy = sd.year
        ey = ed.year
        
        # Months returned as properties of dates are one indexed, (Jan = 1)
        # which complicates the code slightly
        sm = sd.month
        em = ed.month
        
        # find out whether we need to  divide dataset and submit multiple requests
        ds = DATASETS[self.datasetlstbx.curselection()[0]]  # should be a list with a single item
        #print (ds)

        #TODO: refactor the following to avoid code duplication
        #TODO:Refactor to avoid generating a huge data structure

        # I use "in" in the following code because the dispatch method
        # is dependent on a broad class of the data, not the exact identity.
        # This should make it easier to add data types in future.
        if ds in ['GHCND']:
            # requests need to be dispatched by months, to avoid requesting 
            # too much data
            for ayr in range(sy,ey+1):   #Range does not include the last value
                for amo in range(1,13):
                    if (ayr == sy and amo < sm) or (ayr == ey and amo>em):
                        pass
                    else:
                        #We are in the required time period
                        self.themsg.write('Month= %s, year = %s' % (amo, ayr))
                        s,e = getmonthrange(ayr, amo)                    
                        a = self.retreivedata(start = s, end = e, code=ds)
                        if 'results' in a.keys():  # could check for return status
                            if self.datform.get() == "Wide Form":
                                self.appendwidedata(a['results'], self.filename.get())
                            else:
                                self.appendtalldata(a['results'],self.filename.get())
        elif ds in ['GSOM']:
            # Requests need to be dispatched by year to avoid requesting
            # too much data
            for ayr in range(sy,ey+1):   #Range does not include the last value
                self.themsg.write('Year = %s' % (ayr))
                smo = sm if ayr == sy else 1
                emo = em if ayr == ey else 12
                s,e = ('%d-%02d-01' % (ayr,smo), '%d-%02d-31' % (ayr,emo))               
                a = self.retreivedata(start = s, end = e, code=ds)
                if a:
                    if 'results' in a.keys():  # could also check for return status
                        if self.datform.get() == "Wide Form":
                            self.appendwidedata(a['results'], self.filename.get())
                        else:
                            self.appendtalldata(a['results'],self.filename.get())
        elif ds in ['GSOY']:
            # Requests must be submitted for data from less than ten years at a
            # time to avoid error Her we use five years.
            startyr = sy
            while startyr < ey:
                self.themsg.write('Downloading from %d ....' % startyr)
                endyr = startyr+4
                if endyr>ey:
                    endyr=ey
                s,e = ('%d-01-01' % startyr, '%d-12-31' % endyr)
                a = self.retreivedata(start = s, end = e, code=ds, limit=1000)
                if a:
                    if 'results' in a.keys():  # could check for return status
                        if self.datform.get() == "Wide Form":
                            self.appendwidedata(a['results'], self.filename.get())
                        else:
                            self.appendtalldata(a['results'],self.filename.get())
                startyr += 5

#%%
    def appendtalldata(self, dat, savefile = 'talldatafile.csv'):
        '''Save data in "tall" format'''
        openhow = 'w' if self.firstdata else 'a'
        titles = [u'date', u'datatype', u'value', u'attributes'] # set first four data columns
        #TODO  fix this - -assumes data availability is always the same through data
        # which may not always be true
        for item in dat[0]:
            if item in titles:  # Skip any titles we've already added
                continue
            titles.append(item)
    
        with open(savefile, openhow) as thefile:
            csvwriter = csv.DictWriter(thefile, titles,
                                       extrasaction='raise',
                                       lineterminator = '\n')
            if self.firstdata:
                csvwriter.writeheader()
            for item in dat:
                csvwriter.writerow(item)
        self.firstdata=False
                
#%%              
    def appendwidedata(self, dat, outfilename= 'widedatafile.csv'):  
        '''Save Data in "wide" format  
          Assumes data is sorted into blocks by DATE and STATION.
        '''
        openhow = 'w' if self.firstdata else 'a'
        
        #TODO  fix this -- assumes data headers are always the same through data
        # which may not always be true. perhaps make titles a propery of the 
        # downloader? set in the method that dispatches the data download
        # or make  subclass of writer that can modify the dictionary in some
        # sensible way...
        titles = [u'date', u'station'] # set up first data columns
        for item in dat:
            if item['datatype'] in titles:
                continue
            titles += [item['datatype'], item['datatype']+ 'attr']   # Add columns for this data type
           
            
        # actually save the data
        with open(outfilename, openhow) as outfile:
            mydictwriter = csv.DictWriter(outfile, titles,
                                          extrasaction='raise',
                                          lineterminator = '\n')
            if self.firstdata:
                mydictwriter.writeheader()
            thedate = None
            newline = None
            # Every time we hit a line with a new date, we start a new data row
            #Note the assumption here that the data is sorted
            for item in dat:
                if item['date'] != thedate:            # We're starting a new data row.
                    if thedate:                        # If the date exists, then so does the rest of a data row, so record it.
                        mydictwriter.writerow(newline) # First we write the line we were bulding up
                    thedate = item['date']                  # Then we assemble the new line
                    newline = {'date' : thedate[:10]}  # We only want the first 10 characters of the timestamp
                    newline['station'] = self.station.get()
                newline[item['datatype']] = item['value'] # Whether it is a first line for a new date or not, record the data value in the growing line
                if 'attributes' in item.keys():
                    newline[item['datatype']+'attr'] = item['attributes']#Add attributes, Because of commas, these end up quoted
                else:
                    newline[item['datatype']+'attr'] = 'None'
            if thedate:                                # We've fallen off the end of the iterator, so print out the final line -- if there's data ato record
                mydictwriter.writerow(newline)
        self.firstdata = False
#%%
 
    def Quit(self):
        'Quit the Application, after asking to confirm'
        if msgbox.askokcancel("Quit", "Do you really wish to quit?"):
            m = self.master
            self.quit()
            m.destroy()


root = tk.Tk()
#root.minsize(rootMinWindow[0], rootMinWindow[1])
app = NOAAWeatherApp(master=root)
root.protocol("WM_DELETE_WINDOW", app.Quit)

app.mainloop()
