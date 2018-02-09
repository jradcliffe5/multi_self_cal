import re
import sys
import traceback
import logging
import math

def setup_logging_to_file(filename):
    logging.basicConfig( filename='./'+filename,
                         filemode='w',
                         level=logging.DEBUG,
                         format= '%(asctime)s - %(levelname)s - %(message)s',
                       )

def extract_function_name():
    """Extracts failing function name from Traceback
    by Alex Martelli
    http://stackoverflow.com/questions/2380073/\
    how-to-identify-what-function-call-raise-an-exception-in-python
    """
    tb = sys.exc_info()[-1]
    stk = traceback.extract_tb(tb, 1)
    fname = stk[0][3]
    return fname

def log_exception(e):
    logging.error(
    "Function {function_name} raised {exception_class} ({exception_docstring}): {exception_message}".format(
    function_name = extract_function_name(), #this is optional
    exception_class = e.__class__,
    exception_docstring = e.__doc__,
    exception_message = e.message))

def headless(inputfile):
    ''' Parse the list of inputs given in the specified file. (Modified from evn_funcs.py)'''
    INPUTFILE = open(inputfile, "r")
    control = {}
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
            value = value.replace(' ','').strip()
            valuelist = value.split(',')
            if len(valuelist) == 1:
                if valuelist[0] == '0' or valuelist[0]=='1' or valuelist[0]=='2':
                    control[param] = int(valuelist[0])
                else:
                    control[param] = str(valuelist[0])
            else:
                control[param] = ','.join(valuelist)
    return control

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

def dbcon_combine(disk):
	tasav = AIPSTask('TASAV')
	uvavg = AIPSTask('UVAVG')
	avspc = AIPSTask('AVSPC')
	dbcon = AIPSTask('DBCON')
	msort = AIPSTask('MSORT')
	indxr = AIPSTask('INDXR')
	split = AIPSTask('SPLIT')
	splat = AIPSTask('SPLAT')
	uvcop = AIPSTask('UVCOP')
	pca = AIPSCat(indisk)
	print "Your AIPS catalog looks like this:"
	print AIPSCat(indisk)
	print ".... DBCONing....."

	pca = AIPSCat(indisk)
	i = 0
	queueNAME = deque([])
	queueKLAS = deque([])
	queueSEQ  = deque([])
	for fitsfil in pca[indisk]:
		if pca[indisk][i]["klass"] == 'UV' :
			queueNAME.append(pca[indisk][i]["name"])
			queueKLAS.append(pca[indisk][i]["klass"])
			queueSEQ.append(pca[indisk][i]["seq"])
		i = i + 1



	i = 0
	while len(queueNAME)>1 :
		print "Combining " + queueNAME[0] + '.' + queueKLAS[0] + '.' + format(queueSEQ[0]) + " and " + queueNAME[1] + '.' + queueKLAS[1] + '.' + format(queueSEQ[1])
		uvdata1 = AIPSUVData(queueNAME.popleft(), queueKLAS.popleft(), indisk, queueSEQ.popleft())
		uvdata2 = AIPSUVData(queueNAME.popleft(), queueKLAS.popleft(), indisk, queueSEQ.popleft())
		dbcon.indata = uvdata1
		dbcon.in2data = uvdata2
		dbcon.outname = format(i)
		dbcon.outclass = 'DBCON'
		dbcon.outdisk = disk
		dbcon.doarray = 1
		dbcon.fqcenter = -1
		dbcon.go()

		indxr.indata = AIPSUVData(dbcon.outname, dbcon.outclass, indisk, dbcon.outseq)
		indxr.go()

		queueNAME.append(format(i))
		queueKLAS.append('DBCON')
		queueSEQ.append(1)

		i = i + 1
	print "Final sort and index...."
	data = AIPSUVData(queueNAME[0], queueKLAS[0], indisk, queueSEQ[0])
	msort.indata = data
	msort.outdata = data
	msort.go()
	indxr.indata = data
	indxr.go()
	data.rename('COMBO','DBCON',0)
	print "Tidying up your AIPS catalogue."

	pca = AIPSCat(indisk)
	j = 0
	for fitsfil in pca[indisk]:
		if pca[indisk][j]["klass"] == 'FGED' :
			print "Zapping: " + pca[indisk][j]["name"] + '.' + pca[indisk][j]["klass"] + '.' + format(pca[indisk][j]["seq"])
			uvdata1 = AIPSUVData(pca[indisk][j]["name"], pca[indisk][j]["klass"], indisk, pca[indisk][j]["seq"])
			uvdata1.zap()
		if pca[indisk][j]["klass"] == 'DBCON' :
			if pca[indisk][j]["name"].isdigit() :
				print "Zapping: " + pca[indisk][j]["name"] + '.' + pca[indisk][j]["klass"] + '.' + format(pca[indisk][j]["seq"])
				uvdata1 = AIPSUVData(pca[indisk][j]["name"], pca[indisk][j]["klass"], indisk, pca[indisk][j]["seq"])
				uvdata1.zap()
		j = j + 1
