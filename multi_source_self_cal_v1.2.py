# Multi-Source Self Calibration (MSSC)
#    Copyright (C) 2016  Jack Radcliffe
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

####################################################################################
#
# The script uses inputs specified in a file given on the command line.  It loads
# all the data it finds in the same directory and performs MSSC on wide-field VLBI data
#
#
# If the script fails, the first thing to do is run 'source /aips/LOGIN.CSH'.


import os, re, time, datetime, sys, math, fnmatch
from os.path import join, getsize, isfile
from os import listdir
from datetime import date
from collections import deque
import Utilities
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
import math, time, datetime
from numpy import *
import itertools
from time import gmtime, strftime, localtime
ti = time.time()
import numpy as np
import string


def parse_inp(filename):
	''' Parse the list of inputs given in the specified file. (Modified from evn_funcs.py)'''
	INPUTFILE = open(filename, "r")
	control = dict()

	# a few useful regular expressions
	newline = re.compile(r'\n')
	space = re.compile(r'\s')
	char = re.compile(r'\w')
	comment = re.compile(r'#.*')

	# parse the input file assuming '=' is used to separate names from values
	for line in INPUTFILE:
		if char.match(line):
			line = comment.sub(r'', line)
			line = line.replace("'", '')
			(param, value) = line.split('=')

			param = newline.sub(r'', param)
			param = param.strip()
			param = space.sub(r'', param)

			value = newline.sub(r'', value)
			value = value.strip()
			valuelist = value.split(', ')
			control[param] = valuelist

	return control
	


def checkin(control):
	''' convert the control hash to a list of global variables '''

	global userno, indisk, fitsdir, fitsfil, toload, plotdir, fittpdir, msgkill
	global aggr1, max1, aggr2, max2, rho, ncpu, kickoutsigma, filestart, refan, soli, combinIFLLRR, noteles
	global doflag, doimag, imsize, filenam, niter, itercal, fileend, pointcenRA, pointcenDEC, doscal, doload, APhas
	AIPS.userno = int(control.get('userno',[0])[0])
	filenam = int(control.get('filenam',[10])[0])
	doload = int(control.get('doload',[1])[0])
	doscal = int(control.get('doscal',[1])[0])
	indisk = int(control.get('indisk', [0])[0])
	ncpu = float(control.get('ncpu',[8])[0])
	kickoutsigma = float(control.get('kickoutsigma',[1.5])[0])
	doflag = int(control.get('doflag',[1])[0])
	doimag = int(control.get('doimag',[1])[0])
	itercal = int(control.get('self_cal_iter',[3])[0])
	noteles = int(control.get('noteles',[3])[0])
	#Imaging things
	imsize = int(control.get('imsize', [1024])[0])
	refan = int(control.get('ref_ant',[0])[0])
	if imsize:
		imsize = imsize,imsize
	niter = int(control.get('niter', [1024])[0])
	filestart = str(control.get('file_start',[0])[0])
	fileend = str(control.get('file_end',[0])[0])
	pointcenRA = float(control.get('point_cen_RA',[0])[0])
	pointcenDEC = float(control.get('point_cen_DEC',[0])[0])
	soli = int(control.get('solution_interval',[0])[0])
	combinIFLLRR = int(control.get('combineIFLLRR',[0])[0])
	APhas = str(control.get('amplitude_or_phase',[0])[0])
	# Options relating to loading from disk
	fitsdir = control.get('fitsdir', [])
	toload = 0

	# plotdir defaults to the cwd if not specified
	plotdir = control.get('plotdir', [False])[0]
	if (not plotdir):
		plotdir = os.getcwd()
	if (not os.path.isdir(plotdir)):
		print "Error:", plotdir, "does not exist. Check your inputs."
		sys.exit()
	if (not os.access(plotdir, os.W_OK) ):
		print "Error:", plotdir, "is not writable by you. Check your inputs."

	# fittpdir defaults to the cwd if not specified
	fittpdir = control.get('fittpdir', [False])[0]
	if (not fittpdir):
		fittpdir = os.getcwd()
	if (not os.path.isdir(fittpdir)):
		print "Error:", fittpdir, "does not exist. Check your inputs."
		sys.exit()
	if (not os.access(fittpdir, os.W_OK) ):
		print "Error:", fittpdir, "is not writable by you. Check your inputs."

	
	msgkill = int(control.get('msgkill', [-5])[0])

	if (not AIPS.userno) or (not indisk):
		print 'Error: You must set both your AIPS userno and indisk.'
		sys.exit()



# Check for an inputs file
if len(sys.argv)==2 :
	print "Using inputs specified in", sys.argv[1]
	afile = sys.argv[1]
	if os.path.isfile(afile) :
		control = parse_inp(afile)
	else :
		print "Error:" + afile + "does not exist, quitting."
		sys.exit()
else :
	print "Error: no parameter file specified, quitting."
	print "Usage: parseltongue pipeline.py inputs.txt"
	sys.exit()


def findmaxb(uvdata):
	maxbaseline = 0
	antab = uvdata.table('AN',1)
	for row in antab :
		for row2 in antab :
			xsep = row.stabxyz[0] - row2.stabxyz[0]
			ysep = row.stabxyz[1] - row2.stabxyz[1]
			zsep = row.stabxyz[2] - row2.stabxyz[2]
			hypxy =  math.sqrt((xsep * xsep) + (ysep * ysep))
			hypxyz = math.sqrt((zsep * zsep) + (hypxy * hypxy))
			if hypxyz > maxbaseline :
				maxbaseline = hypxyz
	cellsize = (1.22 * (300000000 / uvdata.header.crval[2]) / maxbaseline) / 3.141592 * 180 * 3600 / 5
	print "maxbaseline = ", maxbaseline, "cellsize = ", cellsize
	return cellsize,cellsize

def degreeradecconvert(RA,Dec): # converts degrees ra dec to hh:mm:ss.ss & dd:mm:ss.ss
	inp = ['',RA,Dec]
	del inp[0]

        degRA=float(inp[0])
	degDec=float(inp[1])
	if degRA < 0:
                degRA=degRA+360

        if degRA > 360:
                print `degRA`+": inputs may not exceed 360!\n"
     

        hh=int(degRA/15)
	mm=int((degRA-15*hh)*4)
	ss=(4*degRA-60*hh-mm)*60
        degreeRA = str(string.zfill(`hh`,2)+':'+string.zfill(`mm`,2)+':'+'%10.8f' % ss)

	dd = int(degDec)
	am = int((degDec-dd)*60)
	sec = (((degDec-dd)*60)-am)*60
	degreeDec = str(string.zfill(`dd`,2)+':'+string.zfill(`am`,2)+':'+'%10.8f' % sec)
	return [degreeRA,degreeDec]


def runoffsetradec(uvname): #RA DEC offset from image
			file = open('IMEAN' + uvname[0:8] + '.txt')
			l = []
			for line in file:
				if 'Skypos' in line:
					l.append(line)
			l = l[1]
			ra = l[11:27] #right ascension
			ra = ra.replace(' ',':')
			dec = l[33:-1]
			dec = dec.replace(' ',':')
			return [ra,dec]


def runoffsetpix(uvname): #RA DEC offset from image NEED TO CHANGE IT
			file = open('IMEAN' + uvname + '.txt')
			l = []
			for line in file:
				if 'Maximum=' in line:
					l.append(line)
			l = l[0]
			x = int(l[24:27]) #pixel coords
			y = int(l[29:32])
			return x,y


def maxamplitude(uvname): #RA DEC offset from image NEED TO CHANGE IT
			file = open('IMEAN' + uvname + '.txt')
			l = []
			for line in file:
				if 'Maximum=' in line:
					l.append(line)
			l = l[0]
			x = float(l[9:15])/(10**float(l[18])) #pixel coords

			return x


def radecconvert(a,b,c,d):
    inp=['',a,b,c,d]
    del inp[0]
    if len(inp)==0:
        print" Program to convert an the angular separation between"
        print" two points in the sky"
        print" Type 'angsep.py RA1 Dec1 RA2 Dec2' to calculate the"
        print" angular separation. All coordinates must be of the"
        print" form hh:mm:ss(.ssssssss) or hh mm ss(.ssssssss)"
        print" (Don't mix!).\n"
        sys.exit()
    
    # Find and replace any ":" and "=" from inputs
    newinp=[]
    for x in inp:
        newinp.append(string.replace(x, ":", " "))
    
    inp=newinp
    
    # Find and delete alphanumeric entries like "RA" and "DEC"
    newline=""
    for x in inp:
        newline=newline+" "
        for y in x:
            newline=newline+y
    
    inp=string.split(newline)
    
    newinp=[]
    for x in inp:
        try:
            newinp.append(float(x))
        except ValueError:
            pass
    
    inp=newinp
    
    if len(inp)==4:
        ra1 =string.split(inp[0], ":")
        dec1=string.split(inp[1], ":")
        ra2 =string.split(inp[2], ":")
        dec2=string.split(inp[3], ":")
    elif len(inp)==12:
        ra1 =inp[0:3]
        dec1=inp[3:6]
        ra2 =inp[6:9]
        dec2=inp[9:12]
    else:
        print" Too few or too many parameters."
        sys.exit()
    
    # Calculate angular separation of declinations
    # Convert them into degrees and calulate the difference
    
    # conversion of right ascension 1:
    ra1hh=(float(ra1[0]))*15
    ra1mm=(float(ra1[1])/60)*15
    ra1ss=(float(ra1[2])/3600)*15
    
    ra1deg=ra1hh+ra1mm+ra1ss
    ra1rad=ra1deg*math.pi/180
    
    # conversion of declination 1:
    dec1hh=abs(float(dec1[0]))
    dec1mm=float(dec1[1])/60
    dec1ss=float(dec1[2])/3600
    
    if float(dec1[0]) < 0:
    	dec1deg=-1*(dec1hh+dec1mm+dec1ss)
    else:
    	dec1deg=dec1hh+dec1mm+dec1ss
    
    dec1rad=dec1deg*math.pi/180
    
    # conversion of right ascension 2:
    ra2hh=float(ra2[0])*15
    ra2mm=(float(ra2[1])/60)*15
    ra2ss=(float(ra2[2])/3600)*15
    
    ra2deg=ra2hh+ra2mm+ra2ss
    ra2rad=ra2deg*math.pi/180
    
    # conversion of declination 2:
    dec2hh=abs(float(dec2[0]))
    dec2mm=float(dec2[1])/60
    dec2ss=float(dec2[2])/3600
    
    if float(dec2[0]) < 0:
    	dec2deg=-1*(dec2hh+dec2mm+dec2ss)
    else:
    	dec2deg=dec2hh+dec2mm+dec2ss
    
    dec2rad=dec2deg*math.pi/180
    
    # Delta RA
    deg  = ((ra2rad-ra1rad)*180/math.pi)
    deg_corrected=math.cos(dec1rad)*deg
    hh   =int(deg/15)
    mm   =int((deg-15*hh)*4)
    ss   =(4*deg-60*hh-mm)*60
    
    hh   =int(deg)
    mm   =int((deg-int(deg))*60)
    ss   =((deg-int(deg))*60-mm)*60
    
    # Delta RA corrected for declination (dms format)
    deg_corrected=math.cos(dec1rad)*deg
    hh   =int(deg_corrected)
    mm   =int((deg_corrected-int(deg_corrected))*60)
    ss   =((deg_corrected-int(deg_corrected))*60-mm)*60
    deltara = hh*3600 + mm*60 + ss
    
    # Delta DEC
    deg = ((dec1rad-dec2rad)*180/math.pi)
    hh   =int(deg)
    mm   =int((deg-int(deg))*60)
    ss   =((deg-int(deg))*60-mm)*60
    deltadec = hh*3600 + mm*60 + ss
    
    return [deltara, deltadec]


checkin(control)
pca = AIPSCat(indisk)
AIPSTask.msgkill = msgkill


for file in listdir('./'):
	if file.endswith('tasav.FITS'): ## CHANGE WHEN DATA ARRIVES
		tasavname = str(file)
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "------------------------------------------------------------"
print "              Loading your uv data file                     "
print "------------------------------------------------------------"

if doload == 1:
	AIPSCat().zap()
	uvname = [] #HDF names of files

	for file in listdir('./'):
		if file.startswith(filestart): 
			uvdataname = str(file) #Load the UV data file & image
			#uvdataname = uvdataname[:-1]
			fitld = AIPSTask('FITLD')
			fitld.datain = ('PWD:' + uvdataname) 
			fitld.ncount = 1
			fitld.doconcat = 1
			fitld.clint = 0
			fitld.wtthresh = 0
			fitld.outdisk = indisk
			fitld.digicor = -1
			fitld.go()
		if file.endswith(fileend):
			uvname.append(uvdataname[0:8])
	# makes a list of file prefixes to be used in multi self cal, all data is loaded.



	for i in range(len(uvname)):
			uvdata = AIPSUVData(uvname[i],'SPLAT',indisk,1) #name the uv file in AIPS
			imagedata = AIPSImage(uvname[i],'IIM001',indisk,1)

			nchan = uvdata.header.naxis[2]
			imagr = AIPSTask('IMAGR')
			imagr.nchav = nchan #use imagr to get a clean model!
			imagr.indata = uvdata
			imagr.outname = uvdata.name
			imagr.cellsize[1:] = findmaxb(uvdata)
			imagr.imsize[1:] = imsize
			imagr.nboxes = 1
			imagr.nfield = 1
			imagr.outdisk = indisk
			imagr.uvwtfn = ''
			imagr.niter = niter
			imagr.go()
	
			imagedatacl = AIPSImage(uvname[i],'ICL001',indisk,1)

			imean = AIPSTask('IMEAN')
			imean.indata = imagedatacl
			imean.indisk = indisk
			imean.doprint = 1
			imean.outtext = 'PWD:IMEAN' + uvname[i] + '.txt'
			imean()
		
			uvsub = AIPSTask('UVSUB') #Divide visibilites by clean model!
			uvsub.indata = uvdata
			uvsub.nmaps = 1
			uvsub.in2data = imagedatacl
			uvsub.inver = 1
			uvsub.outdisk = indisk
			uvsub.ncomp[1] = -1000000
			uvsub.opcode = 'DIV'
			uvsub.go()
		
			uvdata = AIPSUVData(uvname[i],'UVSUB',indisk,1)
			#wtmod = AIPSTask('WTMOD') #change weight relative to amplitude adjustments 
			#wtmod.indata = uvdata
			#wtmod.aparm[1] = (maxamplitude(uvname[i])**2)*(10**10)
			#wtmod.outdisk = indisk
			#wtmod.go()
	
			#uvdata = AIPSUVData(uvname[i],'WTMOD',2,1)
			uvdata2 = WizAIPSUVData(uvname[i], 'UVSUB',indisk,1)
			uvdata2.header['crval'][4] = pointcenRA
			uvdata2.header.update()
			uvdata2.header['crval'][5] = pointcenDEC
			uvdata2.header.update()
			uvdata.rename(name='COMBO',klass='UV', seq=0)
	if len(uvname) > 1:
		dbapp = AIPSTask('DBAPP')
		uvdata = AIPSUVData('COMBO','UV',indisk,1)
		dbapp.inname = 'COMBO'
		dbapp.inclass = 'UV'
		dbapp.indisk = indisk
		dbapp.inseq = 2
		dbapp.in2seq = len(uvname)
		dbapp.outdata = uvdata
		dbapp.outdisk = indisk
		dbapp.go()

	uvsrt = AIPSTask('UVSRT')
	uvsrt.sort = 'TB'
	uvsrt.indata = uvdata
	uvsrt.outdata = uvdata
	uvsrt.outdisk = indisk
	uvsrt.go()

	uvdata.zap_table('CL',1)
	indxr = AIPSTask('INDXR')
	indxr.indata = uvdata
	indxr.cparm[1]=360
	indxr.cparm[2] = 360
	indxr.cparm[3] = 0.25
	indxr.go()

	multi = AIPSTask('MULTI') 
	multi.indata = uvdata
	multi.outdisk = indisk
	multi.outname = 'POINT'
	multi.outclass = 'UVDATA'
	multi.go()
	
	uvdata = AIPSUVData('POINT','UVDATA',indisk,1)

	tabed = AIPSTask('TABED')
	tabed.indata = uvdata
	tabed.indisk=indisk
	tabed.outdisk=indisk
	tabed.inext = 'SU'
	tabed.optype = 'REPL'
	tabed.aparm[1:] = 2, 0, 0, 3, 0
	tabed.keystrng = 'POINT'
	tabed.go()


if doscal == 1:
	uvdata = AIPSUVData('POINT','UVDATA',indisk,1)
	calib = AIPSTask('CALIB')
	clcal = AIPSTask('CLCAL')
	snplt = AIPSTask('SNPLT')
	for i in range(1,itercal+1):
		if i >1:
			nchan = uvdata.header.naxis[2]
			imagr = AIPSTask('IMAGR')
			imagr.nchav = nchan #use imagr to get a clean model!
			imagr.indata = uvdata
			imagr.sources[1] = 'POINT'
			imagr.outname = 'POINT'
			imagr.cellsize[1:] = findmaxb(uvdata)
			imagr.imsize[1:] = imsize
			imagr.nboxes = 1
			imagr.nfield = 1
			imagr.outdisk = indisk
			imagr.uvwtfn = ''
			imagr.niter = niter
			imagr.go()
		calib.indata = uvdata
		calib.outdisk = indisk
		calib.gainuse = i
		if i == 1:
			calib.docalib = 0
		if i>1:
			calib.docalib = 2
		if i == 1:
			calib.smodel[1:] = 1, 0, 0, 0
		if i > 1:
			model = AIPSImage('POINT','ICL001',indisk,i-1)
			calib.in2data = model
			calib.calsour[1] = 'POINT'
			calib.in2disk = 2
			calib.nmaps
			calib.ncomp[1] = -1000000
			calib.nmaps = 1
		calib.refant = refan
		calib.solint = soli
		if combinIFLLRR == 1:
			calib.aparm[1:] = 3, 0, 0, 0, 0, 0, 0, 0
		if combinIFLLRR == 2:
			calib.aparm[1:] = 3, 0, 0, 0, 1, 0, 0, 0
		if combinIFLLRR == 3:
			calib.aparm[1:] = 3, 0, 1, 0, 0, 0, 0, 0
		if combinIFLLRR == 4:
			calib.aparm[1:] = 3, 0, 1, 0, 1, 0, 0, 0
		calib.soltype = 'L1'
		calib.solmode = APhas
		calib.snver = i
		calib.go()

		snplt.indata = uvdata
		snplt.indisk = indisk
		snplt.inext = 'SN'
		snplt.pixrange[1:] = -180,180
		snplt.invers = i
		snplt.nplots = 9
		snplt.optype = 'PHAS'
		if combinIFLLRR == 1:
			snplt.opcode = ''
		if combinIFLLRR == 2:
			snplt.opcode = 'ALIF'
		if combinIFLLRR == 3:
			snplt.opcode = 'ALST'
		if combinIFLLRR == 4:
			snplt.opcode = 'ALSI'
		snplt.dotv = -1
		snplt.go()

		clcal.indata = uvdata
		clcal.interpol = 'AMBG'
		clcal.snver = i
		clcal.invers = i
		clcal.gainver = i
		clcal.gainuse = i+1
		clcal.refant = refan
		clcal.go()
	
	lwpla = AIPSTask('LWPLA')
	lwpla.indata = uvdata
	lwpla.indisk = indisk
	lwpla.plver = 1
	if combinIFLLRR == 1:
		lwpla.invers = int(itercal*(math.ceil((noteles*16)/9.0)))
	if combinIFLLRR == 2:
		lwpla.invers = int(itercal*(math.ceil((noteles*2)/9.0)))
	if combinIFLLRR == 3:
		lwpla.invers = int(itercal*(math.ceil((noteles*8)/9.0)))
	if combinIFLLRR == 4:
		lwpla.invers = int(itercal*(math.ceil(noteles/9.0)))
	lwpla.lpen = 1
	lwpla.outfile = 'PWD:Sol_tabl_' + str(itercal) + 'iter_' +str(soli)+ 'sol.ps'
	lwpla.go()

	tasav = AIPSTask('TASAV')
	tasav.indata = uvdata
	tasav.indisk = indisk
	tasav.outname = 'MFSC_SN'
	tasav.outclass = 'TASAV'
	tasav.outdisk = indisk
	tasav.go()

	fittp = AIPSTask('FITTP')
	fittp.indata = AIPSUVData('MFSC_SN','TASAV',indisk,1)
	fittp.indisk = indisk
	fittp.dataout = 'PWD:MFSC_corr_' + str(itercal) + 'iter_' +str(soli)+ 'sol.TASAV'
	fittp.go()
	
	AIPSUVData('MFSC_SN','TASAV',indisk,1).zap()
	for i in range(1,itercal+1):
		uvdata.zap_table('CL',i+1)
		uvdata.zap_table('SN',i)
	if combinIFLLRR == 1:
		for i in range(1,int(itercal*(math.ceil((noteles*16)/9.0)))+1):
			uvdata.zap_table('PL',i)
	if combinIFLLRR == 2:
		for i in range(1,int(itercal*(math.ceil((noteles*2)/9.0)))+1):
			uvdata.zap_table('PL',i)
	if combinIFLLRR == 3:
		for i in range(1,int(itercal*(math.ceil((noteles*8)/9.0)))+1):
			uvdata.zap_table('PL',i)
	if combinIFLLRR == 4:
		for i in range(1,int(itercal*(math.ceil(noteles/9.0)))+1):
			uvdata.zap_table('PL',i)
	
''' # Nugget of gold to shift offsetts to field centres...
uvdata = AIPSUVData(uvname[i],'WTMOD',2,1) #shift data to centre.
centreRADEC = degreeradecconvert(uvdata.header['crval'][4],uvdata.header['crval'][5])
offsetRADEC = runoffsetradec(uvname[i])
offset = radecconvert(centreRADEC[0],centreRADEC[1],offsetRADEC[0],offsetRADEC[1])
uvfix = AIPSTask('UVFIX')
uvfix.indata = uvdata
uvfix.outdisk = indisk
uvfix.shift[1:] = offset[0],offset[1]
uvfix.go()
uvdata = AIPSUVData(uvname[i],'UVFIX',2,1)
print offset, ofsetRADEC, centreRADEC
'''	
		
