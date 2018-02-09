# Multi-Source Self Calibration (MSSC)
#    Copyright (C) 2018  Jack Radcliffe
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
#
# If the script fails, the first thing to do is run 'source /aips/LOGIN.CSH'.


import os, re, time, datetime, sys, math, fnmatch
from os.path import join, getsize, isfile
from datetime import date
import math, time
import itertools
from time import gmtime, strftime, localtime
ti = time.time()
import numpy as np
import string
from MSSC_functions import *
import platform
import logging

### Setup logger
log_name = "%s.log" % os.path.basename(__file__).split('.py')[0]
setup_logging_to_file(log_name)
logging.info('Beginning %s' % os.path.basename(__file__))
print 'Beginning %s' % os.path.basename(__file__)

try:
	from AIPS import AIPS, AIPSDisk
	from AIPSTask import AIPSTask, AIPSList
	from AIPSData import AIPSUVData, AIPSImage, AIPSCat
	from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
	logging.info('ParselTongue found... continuing')
	print 'ParselTongue found... continuing'
except ImportError:
	logging.critical('No Parseltongue found, please re-run with Parseltongue')
	print('No Parseltongue found, please re-run with Parseltongue')
	sys.exit()


# Check for an inputs file
if len(sys.argv)==2 :
	print "Using inputs specified in", sys.argv[1]
	afile = sys.argv[1]
	if os.path.isfile(afile) :
		inputs = headless(afile)
	else :
		print "Error:" + afile + "does not exist, quitting."
		sys.exit()
else :
	print "Error: no parameter file specified, quitting."
	print "Usage: parseltongue multi_source_self_cal.py inputs.txt"
	sys.exit()

### Pull Inputs
### Running options
do_load   = str(inputs['do_load'])
do_self_cal = str(inputs['do_self_cal'])
use_DBAPP = str(inputs['use_DBAPP'])
### File locations
UV_path = str(inputs['UV_path'])
UV_files = str(inputs['UV_files'])
UV_suffix = str(inputs['UV_suffix'])
apply_pbcor = str(inputs['apply_pbcor'])  ## To be used in conjunction with run_pbcor (EVN only)
pbcor_path = str(inputs['pbcor_path'])     ## path to tasav files for PBCOR (EVN only)
do_calib = str(inputs['do_calib'])     ## Apply calibration & split data
## AIPS specific
AIPS.userno = int(inputs['AIPS_userno'])
indisk = int(inputs['AIPS_indisk'])
AIPSTask.msgkill = int(inputs['AIPS_msgkill'])

### MSSC Setup
noteles = int(inputs['noteles'])                  #Number of telescopes
refant = int(inputs['reference_antenna']) #reference antenna used in phase referencing
### UV stacking
imsize = int(inputs['imsize'])   #Image size in pixels.
print imsize
niter = int(inputs['niter'])                    #Clean iterations


pointcenRA, pointcenDEC = str(inputs['pointing_centre']).split(',')


AIPSCat().zap()

if do_load == 'True':
	logging.info('Loading data into AIPS')
	print('Loading data into AIPS')

	if len(UV_files) !=0:
		UV_files = UV_files.split(',')
	else:
		UV_files = []
		for file in os.listdir(UV_path):
			if file.endswith(UV_suffix):
				UV_files = UV_files + [file]
	uvname = [] #HDF names of files

	for i in range(len(UV_files)):
		os.system('rsync -ar --progress %s%s ./' % (UV_path,UV_files[i]))
		fitld = AIPSTask('FITLD')
		fitld.datain = 'PWD:%s' % UV_files[i]
		fitld.ncount = 1
		fitld.doconcat = -1
		fitld.clint = 0
		fitld.wtthresh = 0
		fitld.outdisk = indisk
		fitld.outname = 'A%s' % i
		fitld.outclass = 'LOAD'
		fitld.digicor = -1
		fitld.go()

		uvdata = AIPSUVData('A%s' % i,'LOAD',indisk,1)

		if do_calib == 'True':
			splat = AIPSTask('')

		#imagedata = AIPSImage(uvname[i],'IIM001',indisk,1)

		nchan = uvdata.header.naxis[2]
		imagr = AIPSTask('IMAGR')
		imagr.nchav = nchan #use imagr to get a clean model!
		imagr.indata = uvdata
		imagr.outname = uvdata.name
		imagr.cellsize[1:] = findmaxb(uvdata)
		imagr.imsize[1:] = imsize,imsize
		imagr.uvwtfn='NA'
		imagr.nboxes = 1
		imagr.nfield = 1
		imagr.outdisk = indisk
		imagr.niter = niter
		imagr.go()

		imagedatacl = AIPSImage('A%s' % i,'ICL001',indisk,1)

		imean = AIPSTask('IMEAN')
		imean.indata = imagedatacl
		imean.indisk = indisk
		imean.doprint = 1
		imean.outtext = 'PWD:IMEAN_%s.txt' % UV_files[i]
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

		uvdata = WizAIPSUVData('A%s' % i, 'UVSUB',indisk,1)
		uvdata.header['crval'][4] = pointcenRA
		uvdata.header.update()
		uvdata.header['crval'][5] = pointcenDEC
		uvdata.header.update()
		uvdata.rename(name='COMBO',klass='UV', seq=0)

		os.system('rm %s' % UV_files[i])
	if use_DBAPP == 'True':
		dbapp = AIPSTask('DBAPP')
		uvdata = WizAIPSUVData('COMBO','UV',indisk,1)
		dbapp.inname = 'COMBO'
		dbapp.inclass = 'UV'
		dbapp.indisk = indisk
		dbapp.inseq = 2
		dbapp.in2seq = len(uvname)
		dbapp.outdata = uvdata
		dbapp.outdisk = indisk
		dbapp.go()
		uvdata.rename('COMBO','DBCON',1,1)
	elif use_DBAPP == 'False':
		dbcon_combine(indisk)
		uvdata = WizAIPSUVData('COMBO','DBCON',1,1)
	else:
		print ''

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

if do_self_cal == 'True':
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
			imagr.imsize[1:] = imsize, imsize
			imagr.nboxes = 1
			imagr.docalib = 2
			imagr.nfield = 1
			imagr.outdisk = indisk
			imagr.uvwtfn = 'NA'
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
			calib.in2disk = indisk
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
	tasav.outname = 'MSSC_SN'
	tasav.outclass = 'TASAV'
	tasav.outdisk = indisk
	tasav.go()

	fittp = AIPSTask('FITTP')
	fittp.indata = AIPSUVData('MSSC_SN','TASAV',indisk,1)
	fittp.indisk = indisk
	fittp.dataout = 'PWD:MSSC_corr_' + str(itercal) + 'iter_' +str(soli)+ 'sol.TASAV'
	fittp.go()

	AIPSUVData('MSSC_SN','TASAV',indisk,1).zap()
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
