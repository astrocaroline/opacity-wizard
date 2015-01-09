#!/usr/bin/env python
import scipy as sp
import numpy as np
import pylab as plt
from scipy import interpolate

def getOpacity(filenumber=200, molname='ch4'):
    #  table = np.genfromtxt('opacs/' + molname+'/fort.' + str(int(filenumber)))
    table=np.genfromtxt('http://www.ucolick.org/~cmorley/data/opacities/'+molname+'/fort.' + str(int(filenumber)))
    wl = table[::,0]
    op = table[::,1]
    for i in range(len(op)):#check to see if grid is full of NaNs/zeros. 
        if np.isnan(op[i])  : # replace NaNs with 1e-220 
            op[i] = 1e-220
        if op[i] < 1e-300 or op[i]==0.0:  # replace 0s with 1e-220
            op[i] = 1e-220
    return wl,op

def getPT(log_pressure_layer=1.5, temp_layer=1600):
    # This function takes a pressure and temperature  and finds the bracketing four 	#
    #	PT points in our opacity grid. 
    # Steps: 
    # 	1. find the bracketing temperatures													#
    #	2. find the bracketing pressures													#
    #	3. for value == grid value, use the grid value twice								#
    # Input variables #
    #	pressure_layer : float of pressure, in bar, for atm. layer						 	#
    #	temp_layer : float of temperature in K for atmopshere layer							#
    #	opacityPTlist : list of [index, P[bar], T[K]] for 736/1060 opacity grid				#
    #		should be an increasing (P, T) array. 											#
    # 		program will return an error and quit if the values are outside the table		#
    # Output variables # 
    # list of indices [n1, n2, n3, n4] corresponding to closest 4 points. Index can be 		#
    # repeated if value == grid value. ZERO INDEXED NOT ONE INDEXED							#
    opacityPTlist = np.genfromtxt('opacs/PTgrid1060.txt', skip_header=1) 
    pressure_layer = 10.0**log_pressure_layer
    
    index_opac = opacityPTlist[::,0]
    p_opac = opacityPTlist[::,1]
    t_opac = opacityPTlist[::,2]

    # Find closest temperature first # 		
    for i in range(len(index_opac)):
        if t_opac[i] == temp_layer:  # if the layer temperature is on the opacity grid, use 
            upperT = t_opac[i]		 # that value
            lowerT = t_opac[i]
            break
        elif t_opac[i] > temp_layer:  # find the index where t_opac > temp_layer; that index 
            upperT = t_opac[i]        # is the upper bound, i-1 is the lower bound
            lowerT = t_opac[i-1]
            break
        elif t_opac[len(index_opac)-1] < temp_layer:  # if the layer temperature is higher than our T grid,
            upperT = t_opac[len(index_opac)-1]        # extrapolate from grid edge. index_opac-1 is the upper bound, 
            lowerT = t_opac[len(index_opac)-20] 	  	  # index_opac-2 is the lower bound
            break

    if pressure_layer > np.max(p_opac) : 
    #	print "warning, layer pressure greater than opacity grid, but we're going to extrapolate! "
        for i in range(len(index_opac)):
            if t_opac[i] == upperT: 			
                upperPupperTindex = i  #this should find the last (largest) P point with that T
                lowerPupperTindex = i-1
            elif t_opac[i] == lowerT: #this should find the last (largest) P point with that T
                upperPlowerTindex = i
                lowerPlowerTindex = i-1
        indices = [upperPupperTindex, lowerPupperTindex, upperPlowerTindex, lowerPlowerTindex]


        return indices

    for i in range(len(index_opac)):
        if t_opac[i] == upperT: 			# just for indexes with the right temperature, 
            if p_opac[i] == pressure_layer:	# find  where p_opac is equal to layer pressure
                upperPupperTindex = i
                lowerPupperTindex = i	
                break
            elif p_opac[i] > pressure_layer: #	if no equal point on the grid, find where p_opac
                upperPupperTindex = i		 # becomes > the layer pressure; i-1 is the lower index
                lowerPupperTindex = i-1 	
                break    
    for i in range(len(index_opac)):		# and repeat the above for the lower temperature bound
        if t_opac[i] == lowerT:             
            if p_opac[i] == pressure_layer:
                upperPlowerTindex = i
                lowerPlowerTindex = i 	
                break
            elif p_opac[i] > pressure_layer:
                upperPlowerTindex = i
                lowerPlowerTindex = i-1 	 
                break
    try:  # we tried to catch values off the grid earlier, but due to non-square grid can still get 
              # errors where we're not on the grid. This catches those, returns [0,0,0,0] and a warning. 
        indices = [upperPupperTindex, lowerPupperTindex, upperPlowerTindex, lowerPlowerTindex]
        x= p_opac[indices[0]], p_opac[indices[1]], p_opac[indices[2]], p_opac[indices[3]]
        y= t_opac[indices[0]], t_opac[indices[1]], t_opac[indices[2]], t_opac[indices[3]]
        return indices
    except UnboundLocalError: 
        return [0,0,0,0]



def getMoleculeList(ch4=True,h2o=True,nh3=True,co=True,co2=True,feh=False,h2s=True,ph3=True,tio=False,vo=False, **args):
	#print args
	molList = []
	for i in args: 
	#	print i
	#	print args[i]
		if i == u'H\u2082O' and args[i]==True: 
			molList.append('h2o')
		if i == u'CH\u2084' and args[i]==True: 
			molList.append('ch4')
		if i == u'CO' and args[i]==True: 
			molList.append('co')
		if i == u'CO\u2082' and args[i]==True: 
			molList.append('co2')
		if i == u'NH\u2083' and args[i]==True: 
			molList.append('nh3')
		if i == u'PH\u2083' and args[i]==True: 
			molList.append('ph3')
		if i == u'FeH' and args[i]==True: 
			molList.append('feh')
		if i == u'TiO' and args[i]==True: 
			molList.append('tio')
		if i == u'VO' and args[i]==True: 
			molList.append('vo')
		if i == u'H\u2082S' and args[i]==True: 
			molList.append('h2s')
	return molList

def getMetallicity(metallicity):
	return metallicity

def interpolateOpacityFiles(pressure_layer, temp_layer, indices, opacity_array):
	# This function takes a set of four indices which represent the pressure and 			#
	# temperatures to interpolate, the true PT point we want, and absorber name and returns	#
	# the array of freq vs. opacity for that absorber. 										#
	# Steps: 
	# 	1. Find P and T of each index in indices											#
	#	2. Interpolate linearly in log10(pressure)											#
	#	3. Then interpolate those linearly in temperature									#
	# (this is an algorithm for a standard bilinear interpolation, 							#
	#	e.g. http://en.wikipedia.org/wiki/Bilinear_interpolation )							#
	# Input variables #
	#	log_pressure_layer : float, the pressure of the layer we're interpolating for in bar	#
	#	temp_layer : float, the temperature of the layer we're intepolating for in K		#
	#	opacityPTlist : list of [index, P[bar], T[K]] for 736/1060 opacity grid				#
	#	indices : ZERO INDEXED indices that we want. Note that opacityPTlist[i][0] will		#
	# 			return i+1. Pressure[i] = opacityPTlist[i][1]								#
	# 	 		[upperPupperTindex, lowerPupperTindex, upperPlowerTindex, lowerPlowerTindex]#
	#	opacity_array : numpy array of opacities for each PT point in our 4 indexed points	#
	# 					length is 4 x len(freqs). 
	# Output variables # 
	# 	opacities_interp : 	the interpolated opacities. len(freqs) 							#
	
	#Get the pressures and the temperatures of our indexed points
    opacityPTlist = np.genfromtxt('opacs/PTgrid1060.txt', skip_header=1) 
    
   
    p1 = opacityPTlist[indices[1]][1]
    p2 = opacityPTlist[indices[0]][1]
    t1 = opacityPTlist[indices[2]][2]
    t2 = opacityPTlist[indices[0]][2]
#	print 'pressures:' , p1,p2
#	print 'temps: ', t1,t2
#	print 'point wanted: ', pressure_layer, temp_layer
	#Make sure our grid is rectangular. (x1,y1), (x2,y1), (x1,y2), (x2,y2)
    if opacityPTlist[indices[3]][1] != p1: sys.exit("bilinear interpolation doesn't have a square grid in interpolateOpacityFiles!")
    if opacityPTlist[indices[2]][1] != p2: sys.exit("bilinear interpolation doesn't have a square grid in interpolateOpacityFiles!")
    if opacityPTlist[indices[3]][2] != t1: sys.exit("bilinear interpolation doesn't have a square grid in interpolateOpacityFiles!")
    if opacityPTlist[indices[1]][2] != t2: sys.exit("bilinear interpolation doesn't have a square grid in interpolateOpacityFiles!")
    
    p1 = np.log10(p1)
    p2 = np.log10(p2)
    t1 = np.log10(t1)
    t2 = np.log10(t2)
    opacity_array = np.log10(opacity_array)
    temp_layer = np.log10(temp_layer)
    pressure_layer = np.log10(pressure_layer)
    #This line fixes a problem where we were extrapolating negative infinity opacities: 
    opacity_array[opacity_array < -500] = -500
    
    #If p1=p2 and t1=t2, no interpolation, return point
    if ((p1==p2) and (t1==t2)):
        return 10.**opacity_array[0]
    #If p1=p2, do linear interpolation in t 
    elif (p1==p2) :
        opacities_interp =  opacity_array[2]  + (opacity_array[0] -opacity_array[2]) * (temp_layer-t1)/(t2-t1)
        return 10.**opacities_interp
	#If t1=t2, do linear interpolation in p 
    elif (t1==t2) :
        opacities_interp =  opacity_array[1]  + (opacity_array[0] -opacity_array[1]) * (pressure_layer-p1)/(p2-p1)
        return 10.**opacities_interp
	
	#Otherwise, do the full bilinear interpolation
#	print 'Q11', opacity_array[3] 
#	print 'Q12', opacity_array[2] 
#	print 'Q21', opacity_array[1] 
#	print 'Q22', opacity_array[0] 
#	print 'calculating R1: ', p2, pressure_layer, p1, opacity_array[3], opacity_array[2]
    R1 = ((p2-pressure_layer) / (p2-p1)) * opacity_array[3] + ((pressure_layer-p1)/(p2-p1)) * opacity_array[2]
    R2 = ((p2-pressure_layer) / (p2-p1)) * opacity_array[1] + ((pressure_layer-p1)/(p2-p1)) * opacity_array[0]
#	print 'R1', R1
#	print 'R2', R2
    opacities_interp = ((t2-temp_layer) / (t2-t1)) * R1 + ((temp_layer-t1)/(t2-t1)) * R2
#	print 'opacities_interp', opacities_interp
    opacities_interp = 10.0**opacities_interp	
    return opacities_interp

def getAbundances(pressures, temps, moleculelist, metallicity='0.0'):
    # This function takes pressure and temperature arrays for the model layers and a list   #
    # of additional molecules that we'd like to include the abundances of. The equilibrium  #
    # abundances of the additional molecules will be returned in the same order as the list.#
    # Input variables #
    #   pressures : numpy array of the pressures of model layers in bar, length N_layer     #
    #   temps : numpy array of the temperatures of model layers in K, length N_layer        #
    #   moleculelist : list of strings of additional molecules we're including              #
    # Output variables # 
    #   molecule_abunds : list of numpy arrays of additional molecules.                     #
    #                        [len(moleculelist) x N_layer]                                  #
    
    #read in the abunds file and the corresponding PT grid
    abundsfile = 'abunds/abunds.'+metallicity
    ptfile = 'abunds/PTpoints.txt'
    f= open(abundsfile, 'r') #opens the file to read in
    firstline = f.readline().split() #reads and splits first line
    nMolecules = len(firstline) #finds the number of "molecules" (also includes atoms etc.)
    abundsGrid = np.genfromtxt(abundsfile, skip_header=1) #get grid of abundances
    f.close()
 
    
    f= open(ptfile, 'r') #opens the file to read in
    firstline = f.readline().split() #reads and splits first line
    PTgrid = np.genfromtxt(ptfile, skip_header=1) #get grid of PTs for abundances
        # the number of lines here corresponds to the number of lines in the abunds file
        # (736 or 1060 usually). These are just written in order here. So the first line
        # in the abundances file corresponds to the first line here in terms of pressure
        # and temperature.
        # Grid is in 3 columns. First is line (file) number. Second pressure in bar. 
        # Third is temperature in K. 
    nlines_ptgrid = len(PTgrid)
    f.close()
    
    pressures = np.array([pressures])
    temps = np.array([temps])
    
    if (nlines_ptgrid != len(abundsGrid)): 
        print "ERROR! Abundance file has diff. line length than PTpoints file!"
    
    #Clean up arrays, replace NaNs with 1e-220 and zeros with 1e-220    
    pressArr = PTgrid[:,1] #pressures array, same length as abunds file 
    tempArr = PTgrid[:,2]  #temperature array, same length as abunds file    
    for i in range(len(abundsGrid)):#check to see if grid is full of NaNs. 
        for j in range(len(abundsGrid[i])):
            if np.isnan(abundsGrid[i][j])  : # replace NaNs with 1e-220 
                abundsGrid[i][j] = 1e-220
            if abundsGrid[i][j] == 0.0:  # replace 0s with 1e-220
                abundsGrid[i][j] = 1e-220
    
    #dictionary of columns          
    columns = {'e-': 2, 'h2':3, 'h':4, 'h+':5,'h-':6 , 'h2-':7, 'h2+':8, 'h3+':9, 'he':10,
                'h2o':11, 'ch4':12, 'co':13, 'nh3':14, 'n2':15, 'ph3':16,  'h2s':17, 'tio':18,
                'vo':19, 'fe':20, 'feh':21, 'crh':22, 'na':23, 'k':24, 'rb':25, 'cs':26, 'co2':27}
    
    #setup array for abundances
    molecule_abunds = np.zeros((len(moleculelist), len(pressures)))
    
    for i in range(len(moleculelist)):
        moleculename = moleculelist[i]
        try: 
            columnnum = columns[moleculename]
        except KeyError:
            print 'molecule name: ', moleculename, 'not found in abunds dictionary!'
 
        abunds_molecule_grid = abundsGrid[::,columnnum]
        
        #set up meshes of our variables
        mT,mP = np.meshgrid(temps, np.log10(pressures))
        #print mT
        mT= np.transpose(mT)
        mP= np.transpose(mP) # meshgrid creates a grid of all the P, T points. We then tranpose
                            # them to be the way around that we prefer. 
        #print mT
    
        molecule_abunds[i] = 10.0**np.diag(interpolate.griddata((tempArr, np.log10(pressArr)), 
                        np.log10(abunds_molecule_grid), (mT, mP), method='linear'))
                        # interpolate the abunds grid onto the model P, T grid made earlier. 
                        # Then takes just the diagonal for just the P, T points we care about, 
                        # e.g. (1,1) gives us the top point in the atmosphere. 
 
        for j in range(len(molecule_abunds[i])):
            if np.isnan(molecule_abunds[i][j]): 
                #   print i, j, molecule_abunds[i][j]
                    #fill in NaN missing numbers with previous layer's value
                molecule_abunds[i][j]=molecule_abunds[i][j-1]
 
    return molecule_abunds 
    
def makeOpacityPlot(opac_kwargs, opac_result, abund_result, met_result):
    wl,op=getOpacity(filenumber=1)
    op_array = np.zeros((4,len(op)))
    #color_array=['Blue', 'Red', 'DarkViolet', 'LimeGreen', 'DarkOrange', 'DeepPink', 'Teal']
    color_list = plt.cm.rainbow(np.linspace(0, 1, len(abund_result) ))
    leg_dict = { 'h2o':'H$_2$O', 'ch4':'CH$_4$', 'co':'CO', 'nh3':'NH$_3$', 'n2':'N$_2$', 'ph3':'PH$_3$',  'h2s':'H$_2$S', 'tio':'TiO',
                'vo':'VO', 'feh':'FeH', 'crh':'CrH', 'na':'Na', 'k':'K', 'rb':'Rb', 'cs':'Cs', 'co2':'CO$_2$'}

    fig1 = plt.figure(num=1, figsize=(12,8))

    for j in range(len(abund_result)):
     #   print abund_result[j]
        for i in range(len(opac_result)):
            ii = opac_result[i]+1
            wl,op_array[i] = getOpacity(filenumber=ii, molname=abund_result[j])
        op = interpolateOpacityFiles(10.0**opac_kwargs['log_pressure_layer'], opac_kwargs['temp_layer'], opac_result, op_array)
        abunds = getAbundances(10.0**opac_kwargs['log_pressure_layer'], opac_kwargs['temp_layer'], [abund_result[j]], metallicity=met_result)
        plt.loglog(wl,op*abunds[0],color=color_list[j],label=leg_dict[abund_result[j]])

    plt.xlim(0.8,20);
    plt.ylim(1e-35,1e-18);
    plt.xlabel('wavelength ($\mu$m)', size='xx-large');
    plt.ylabel('opacity $\\times$ mixing ratio (cm$^2$/molecule) ', size='xx-large');
    plt.minorticks_on();
    plt.tick_params(length=10, width=1, labelsize='x-large', which='major');
    plt.tick_params(length=5, width=1,  which='minor');
    plt.legend(loc='upper left',frameon=False, ncol=len(abund_result)/3);
    plt.annotate('[M/H]: ' + met_result, xy=(0.7,0.92), xycoords='axes fraction', size='x-large', color='Black')
    plt.annotate('pressure: ' + str(round(10.0**opac_kwargs['log_pressure_layer'],1)) + ' bar', xy=(0.7,0.87), xycoords='axes fraction', size='x-large', color='Black')
    plt.annotate('temperature: ' + str(opac_kwargs['temp_layer']) + ' K', xy=(0.7,0.82), xycoords='axes fraction', size='x-large', color='Black')
	

    return fig1

def makeAbundsPlot(opac_kwargs, opac_result, abund_result, met_result): 
	abunds = getAbundances(10.0**opac_kwargs['log_pressure_layer'], opac_kwargs['temp_layer'], abund_result, metallicity=met_result)
	color_list = plt.cm.rainbow(np.linspace(0, 1, len(abund_result) ))
	leg_dict = { 'h2o':'H$_2$O', 'ch4':'CH$_4$', 'co':'CO', 'nh3':'NH$_3$', 'n2':'N$_2$', 'ph3':'PH$_3$',  'h2s':'H$_2$S', 'tio':'TiO',
                'vo':'VO', 'feh':'FeH', 'crh':'CrH', 'na':'Na', 'k':'K', 'rb':'Rb', 'cs':'Cs', 'co2':'CO$_2$'}
	xaxis = np.arange(0,len(abunds),1)
	molnames=[]
	fig = plt.figure(num=1, figsize=(11,7))
	ax = plt.subplot(111)
	for i in range(len(abunds)):
		rect1=ax.bar(0.5+i, abunds[i], width=0.8, color=color_list[i],edgecolor=color_list[i])
		molnames.append(leg_dict[abund_result[i]])
		rect=rect1[0]
		height=rect.get_height()
		plt.ylim(0,max(abunds)+0.15*max(abunds))
		ax.text(rect.get_x()+rect.get_width()/2., height+0.025*max(abunds), '%.2e'%height, ha='center', va='bottom', color=color_list[i], size='large')
	spines_to_remove = ['top', 'right']
	for spine in spines_to_remove:
		ax.spines[spine].set_visible(False)
	ax.xaxis.set_ticks_position('none')
	ax.yaxis.set_ticks_position('none')
	ax.yaxis.grid(color='white', linestyle='solid', lw=1.5)
	plt.xticks(xaxis+0.9, molnames, size='xx-large')
	ax.tick_params(axis='y', which='major', labelsize='x-large')
	plt.ylabel('number mixing ratio', size='xx-large')	
	plt.annotate('[M/H]: ' + met_result, xy=(0.7,0.95), xycoords='axes fraction', size='xx-large', color='Black')
	plt.annotate('pressure: ' + str(round(10.0**opac_kwargs['log_pressure_layer'],1)) + ' bar', xy=(0.7,0.9), xycoords='axes fraction', size='xx-large', color='Black')
	plt.annotate('temperature: ' + str(opac_kwargs['temp_layer']) + ' K', xy=(0.7,0.85), xycoords='axes fraction', size='xx-large', color='Black')
	
	return fig

