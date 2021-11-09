
# Module de dessin d'un karyotype

import sys
import itertools

from . import myTools
from . import myPsOutput

myTools.addModuleOptions("karyo", [("roundedChr",bool,False), ("resolution",int,1), ("showText",bool,True), ("drawBorder",bool,False), ("defaultColor",str,""), ("penColor",str,"black")] )


def drawKaryo(data, arguments, x0=1, y0=1, lx=None, ly=None, dx=None, dy=None, bysize=False):

	defaultColor = arguments["karyo:defaultColor"] if len(arguments["karyo:defaultColor"]) > 0 else None

	# Gestion de plusieurs pistes par chromosome
	data = [(c,l if isinstance(l, tuple) else [l]) for (c,l) in data]
	if bysize:
		data.sort(key=lambda x: len(x[1][0]), reverse=True)

	# Les dimensions doivent etre specifiees
	assert not ((lx is None) and (dx is None))
	assert not ((ly is None) and (dy is None))
	
	if arguments["karyo:showText"]:
		if ly is not None:
			ly -= 1
		y0 += 1
	
	if dx is None:
		dx = (lx*3.)/(5.*len(data)-2.)
	
	if dy is None:
		dy = ly / float(max([len(x[0]) for (_,x) in data]))

	# On dessine
	for (c,all) in data:
		size = len(all[0])
		newall = []
		for currl in all:
			newcurr = []
			for (col,items) in itertools.groupby(currl):
				ly = len(list(items))
				if ly >= arguments["karyo:resolution"]:
					newcurr.extend( [col]*ly )
			newall.append(newcurr)
		all = newall
		size = max(len(x) for x in all)

		def printBorder():
			print("newpath")
			if arguments["karyo:roundedChr"]:
				print("%.5f cm %.5f cm %.5f cm 180 360 arc" % (x0+dx/2., y0+dx/2., dx/2.))
				print("%.5f %.5f 2cm rlineto" % (0,size*dy-dx))
				print("%.5f cm %.5f cm %.5f cm 0 180 arc" % (x0+dx/2., y0+size*dy-dx/2., dx/2.))
				print("%.5f %.5f 2cm rlineto" % (0,-size*dy+dx))
			else:
				print("%.5f %.5f 2cm moveto" % (x0,y0))
				print("%.5f %.5f 2cm rlineto" % (0,size*dy))
				print("%.5f %.5f 2cm rlineto" % (dx,0))
				print("%.5f %.5f 2cm rlineto" % (0,-size*dy))
			print("closepath")

		if arguments["karyo:roundedChr"]:
			printBorder()
			print("clip")

		print("0 cm setlinewidth")
		nbl = float(len(all))
		for (i,currl) in enumerate(all):
			y = y0
			for (col,items) in itertools.groupby(currl):
				ly = len(list(items))
				if ly < arguments["karyo:resolution"]:
					continue
				if col is None:
					col = defaultColor
				ly *= dy
				myPsOutput.drawBox(x0+i*dx/nbl, y, dx/nbl, ly, col, col)
				y += ly
		print("0.01 cm setlinewidth")
		
		if arguments["karyo:drawBorder"]:
			myPsOutput.setColor(arguments["karyo:penColor"], "color")
			printBorder()
			for i in range(1,len(all)):
				print("%.5f %.5f 2cm moveto" % (x0+i*dx/nbl,y0))
				print("%.5f %.5f 2cm rlineto" % (0,size*dy))
			print("stroke")
		
		if arguments["karyo:roundedChr"]:
			print("initclip")

		if arguments["karyo:showText"] and (c is not None):
			myPsOutput.drawText(x0, y0-1, str(c), arguments["karyo:penColor"])
		
		x0 += (5.*dx)/3.


