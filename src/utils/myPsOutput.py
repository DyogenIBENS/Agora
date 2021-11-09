
# Module d'ecriture dans un fichier PostScript
import os

# Les petits noms que je donne a mes couleurs (nom -> nom UNIX)
colorTransl = {}
# Les definitions des couleurs (valeurs RGB -> nom UNIX) pour eviter de creer plusieurs fois la meme couleur
colorTableRGB2UNIX = {}
# La table inverse
colorTableUNIX2RGB = {}


#
# L'en-tete PostScript
#
def printPsHeader(landscape = False):
	print(strPsHeader)
	print(strMyDef)
	print(strFontDef)
	print("1 setlinejoin")
	print("0.001 cm setlinewidth")
	print()

	initColor()

	if landscape:
		print("90 rotate")
		print("0 -21 cm translate")
		print()
		return (29.7,21.)

	print()
	return (21.,29.7)


#
# Le pied de page PostScript
#
def printPsFooter():
	print("showpage")



#
# Charge le fichier de definitions des couleurs en RGB
# Initialise toutes les couleurs couramment utilisees
#
def initColor(silent = False):

	knownLocations = [
		"/etc/X11/rgb.txt",
		"/opt/X11/share/X11/rgb.txt",
	]
	location = [f for f in knownLocations if os.path.exists(f)][0]
	# La liste des couleurs et leurs valeurs RGB
	f = open(location, 'r')
	for l in f:
		if l.startswith("!"):
			continue
		c = l.split()
		s = "".join(c[3:]).lower()
		if s in colorTableUNIX2RGB:
			continue
		rgb = tuple(int(x) for x in c[:3])
		if not silent:
			print("/%s [%.3f %.3f %.3f] def" % (s, rgb[0]/255., rgb[1]/255., rgb[2]/255.))
		colorTableRGB2UNIX[rgb] = s
		colorTableUNIX2RGB[s] = rgb

	f.close()


	c5 = []
	c5 += ["red4", "coral4", "firebrick", "red2", "coral2", "darkorange", "gold"]
	c5 += ["yellow", "khaki1","wheat1","peachpuff","lightsalmon", "hotpink2", "magenta2", "darkorchid2", "purple2", "darkorchid4"]
	c5 += ["blue2", "royalblue2", "blue4"]
	c5 += ["turquoise4", "darkseagreen4", "chartreuse4","mediumaquamarine",(121,204,61), "chartreuse2", "olivedrab2", "darkolivegreen1"]
	#c5 += ["turquoise4", "darkseagreen4", "chartreuse4","mediumaquamarine", "chartreuse2", "olivedrab2", "darkolivegreen1"]
	c5 += ["darkseagreen1", "paleturquoise2", "lightblue", "skyblue1", "turquoise2", "lavender","thistle2"]
	c5 += [(204,204,153),"lightgoldenrod3","ivory2","honeydew3","slategray"]
	#c5 += ["lightgoldenrod3","ivory2","honeydew3","slategray"]

	lightColors = ["salmon", "PaleTurquoise2", "DarkSeaGreen1","khaki1", "thistle2", "PeachPuff","LightBlue", "SkyBlue1","LightGoldenrod3", "wheat1",  "DarkOliveGreen1", "lavender"]

	darkColors = ["red1", "turquoise2", "DarkGreen", "yellow", "coral2", "OliveDrab2", "orange", "MediumAquamarine", "blue2", "firebrick4", "LightSalmon", "DarkViolet", "magenta2", "DarkSeaGreen4", "DarkSlateBlue", "yellow4", "grey62", "gold", "PeachPuff2", "HotPink4", "firebrick", "purple4"]
	darkColors = ["red3", "chartreuse4", "turquoise2", "gold", "coral2", "OliveDrab2", "orange", "DarkSeaGreen4", "firebrick4", "RoyalBlue4", "LightSalmon", "DarkViolet", "magenta2", "MediumAquamarine"]
	#darkColors = ["turquoise2", "gray85", "yellow", "coral2", "OliveDrab2", "MediumAquamarine", "blue2", "LightSalmon", "magenta2", "DarkSeaGreen2", "PaleTurquoise3", "khaki2", "thistle3", "PeachPuff2", "HotPink2", "firebrick", "purple2"]

	greekLetters = ["ALPHA", "BETA", "DELTA", "EPSILON", "GAMMA", "PHI"]
	craniateColors = ["DarkOrange", "RoyalBlue4", "chartreuse4", "gold", "DarkOrchid4", "red3"]
	#craniateColors = ["DarkOrange", "RoyalBlue2", "chartreuse2", "gold", "DarkOrchid1", "red2"]

	# Les couleurs claires sont pour les nombres negatifs et les lettres minuscules
	ordre = lightColors + craniateColors + darkColors
	for (i,c) in enumerate(ordre):
		colorTransl[str(-(i+1))] = c.lower()
		if i < 26:
			colorTransl[chr(i+97)] = c.lower()

	# Les couleurs foncees sont pour les nombres positifs et les lettres majuscules
	ordre = darkColors + lightColors + craniateColors
	#ordre = craniateColors + darkColors + lightColors
	for (i,c) in enumerate(ordre):
		colorTransl[str(i+1)] = c.lower()
		if i < 26:
			colorTransl[chr(i+65)] = c.lower()

	for (i,c) in enumerate(greekLetters):
		colorTransl[c] = craniateColors[i]


########################
# Gestion des couleurs #
########################

def getLinearGradient(colors, nelem):
	l = []
	nc = (nelem-1.) / (len(colors)-1.)
	for i in range(nelem-1):
		ip = int(i/nc)
		l.append( alphaColor(colors[ip], colors[ip+1], i/nc-ip) )
	l.append(colors[-1])
	return l


def getCubicGradient(colors, nelem):
	import myMaths
	l = len(colors)
	tmp = (l-1.)/(nelem-1.)
	interpol = myMaths.myInterpolator.getMultDim(myMaths.myInterpolator.oneDimCubic, list(range(l)), colors)
	return [interpol(i*tmp) for i in range(nelem)]


def alphaColor(xxx_todo_changeme, xxx_todo_changeme1, alpha):
	(r1,g1,b1) = xxx_todo_changeme
	(r2,g2,b2) = xxx_todo_changeme1
	return (int(round(r1*(1.-alpha)+r2*alpha)), int(round(g1*(1.-alpha)+g2*alpha)), int(round(b1*(1.-alpha)+b2*alpha)))

def YUV2RGB(xxx_todo_changeme2):
	(Y,U,V) = xxx_todo_changeme2
	R = Y + 1.140*V
	G = Y - 0.395*U - 0.581*V
	B = Y + 2.032*U

	R = int(R*255)
	G = int(G*255)
	B = int(B*255)

	return (R,G,B)



#######################
# Fonctions de dessin #
#######################


# Les coordonnes sont en cm (portrait/paysage)
#   X: de gauche (0) a droite (21/29.7)
#   Y: de bas (0) en haut (29.7/21)
# Les couleurs sont facultatives, utiliser None pour ne pas en specifier.
#   Sinon, utiliser le nom UNIX (cf getColor pour avoir ce nom)
# Les angles sont en degres

def setColor(C, txt):

	# Un nom comme on les aime
	if C in colorTableUNIX2RGB:
		s = C
	# Des raccourcis
	elif str(C) in colorTransl:
		s = colorTransl[str(C)]
	# Un triplet (r,g,b)
	else:
		try:
			if len(C) == 3:
				(r,g,b) = list(map(int,C))
			elif C[0] == '#':
				try:
					(r,g,b) = [int(x) for x in C[1:].split(':')]
				except ValueError:
					(r,g,b) = [int(x,16) for x in C[1:].split(':')]
			else:
				return
		except TypeError:
			return

		if (r,g,b) in colorTableRGB2UNIX:
			s = colorTableRGB2UNIX[(r,g,b)]
		else:
			s = "tmp%d" % len(colorTableUNIX2RGB)
			colorTableRGB2UNIX[(r,g,b)] = s
			colorTableUNIX2RGB[s] = (r,g,b)
			print("/%s [%.3f %.3f %.3f] def" % (s, float(r)/255., float(g)/255., float(b)/255.))
	print(s, txt)


def drawLine(X, Y, L, H, C):
	setColor(C, "color")
	print("%.6g %.6g %.6g %.6g myline" % (L,H, X,Y))


def drawBox(X, Y, L, H, Cb, Cr):
	print("%.6g %.6g %.6g %.6g mybox" % (L,H, X,Y))
	setColor(Cr, "myfill")
	setColor(Cb, "color stroke")


def drawCross(X, Y, L, H, C):
	print("newpath")
	setColor(C, "color")
	print("%.6g cm %.6g cm moveto" % (X, Y))
	print("%.6g cm %.6g cm rlineto" % (L, H))
	print("%.6g cm %.6g cm moveto" % (X+L, Y))
	print("%.6g cm %.6g cm rlineto" % (-L, H))
	print("closepath")
	print("stroke")

def drawCircle(X, Y, R, A, B, Cb, Cr):
	print("newpath")
	setColor(Cb, "color")
	print("%.6g cm %.6g cm %.6g cm %.6g %.6g arc" % (X, Y, R, A, B))
	setColor(Cr, "myfill")
	print("stroke")

def drawArrowR(X, Y, L, H, P, Cb, Cr):
	print("newpath")
	setColor(Cb, "color")
	print("%.6g cm %.6g cm moveto" % (X, Y))
	print("%.6g cm 0 cm rlineto" % L)
	print("%.6g cm %.6g cm rlineto" % (P, H/2))
	print("%.6g cm %.6g cm rlineto" % (-P, H/2))
	print("%.6g cm 0 cm rlineto" % (-L))
	print("closepath")
	setColor(Cr, "myfill")
	print("stroke")

def drawArrowL(X, Y, L, H, P, Cb, Cr):
	print("newpath");
	setColor(Cb, "color")
	print("%.6g cm %.6g cm moveto" % (X, Y+(H/2)))
	print("%.6g cm %.6g cm rlineto" % (P, H/2))
	print("%.6g cm 0 cm rlineto" % L)
	print("0 cm %.6g cm rlineto" % (-H))
	print("%.6g cm 0 cm rlineto" % (-L))
	print("closepath")
	setColor(Cr, "myfill")
	print("stroke")

def drawText(X, Y, T, C):
	setColor(C, "color")
	print("(%s) %.6g %.6g mytext" % (T, X,Y))



strPsHeader = """%!PS-Adobe-3.0
%%DocumentData: Clean7bit
%%Creator: myPsOutput
%%PageOrder: Ascend
%%Pages: 1
%%DocumentFonts: Helvetica
%%LanguageLevel: 1
%%EndComments"""

strMyDef = """
/color { aload pop setrgbcolor } def
/cm {28.3464567 mul} def
/2cm { cm exch cm exch } def
/myline { newpath 2cm moveto 2cm rlineto closepath stroke } def
/mybox { newpath 2cm moveto 2cm exch dup 0 rlineto exch 0 exch rlineto neg 0 rlineto closepath } def
/myfill { gsave color fill grestore } def
/mytext { 2cm moveto show } def
"""

strFontDef = """
%%BeginProlog

% Font encoding-vector for iso-8859-1 (latin-1)
/font_encoding_vector [
/.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef 
/.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef 
/.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef 
/.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef 
/space /exclam /quotedbl /numbersign /dollar /percent /ampersand /quoteright 
/parenleft /parenright /asterisk /plus /comma /hyphen /period /slash 
/zero /one /two /three /four /five /six /seven /eight /nine /colon 
/semicolon /less /equal /greater /question /at /A /B /C /D /E /F /G 
/H /I /J /K /L /M /N /O /P /Q /R /S /T /U /V /W /X /Y /Z /bracketleft 
/backslash /bracketright /asciicircum /underscore /quoteleft /a /b /c 
/d /e /f /g /h /i /j /k /l /m /n /o /p /q /r /s /t /u /v /w /x /y /z 
/braceleft /bar /braceright /asciitilde /.notdef 
/.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef 
/.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef 
/.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef 
/.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef 
/space /exclamdown /cent /sterling /currency /yen /brokenbar /section 
/dieresis /copyright /ordfeminine /guillemotleft /logicalnot /hyphen 
/registered /macron /degree /plusminus /twosuperior /threesuperior /acute 
/mu /paragraph /bullet /cedilla /dotlessi /ordmasculine /guillemotright 
/onequarter /onehalf /threequarters /questiondown /Agrave /Aacute 
/Acircumflex /Atilde /Adieresis /Aring /AE /Ccedilla /Egrave /Eacute 
/Ecircumflex /Edieresis /Igrave /Iacute /Icircumflex /Idieresis /Eth 
/Ntilde /Ograve /Oacute /Ocircumflex /Otilde /Odieresis /multiply /Oslash 
/Ugrave /Uacute /Ucircumflex /Udieresis /Yacute /Thorn /germandbls 
/agrave /aacute /acircumflex /atilde /adieresis /aring /ae /ccedilla 
/egrave /eacute /ecircumflex /edieresis /igrave /iacute /icircumflex 
/idieresis /eth /ntilde /ograve /oacute /ocircumflex /otilde /odieresis 
/divide /oslash /ugrave /uacute /ucircumflex /udieresis /yacute /thorn 
/ydieresis 
] def

/MF {   % fontname newfontname -> -     make a new encoded font
  /newfontname exch def
  /fontname exch def
  /fontdict fontname findfont def
  /newfont fontdict maxlength dict def
  fontdict { exch
    dup /FID eq { pop pop } 
      {  % copy to the new font dictionary
         exch newfont 3 1 roll put } ifelse
  } forall
  newfont /FontName newfontname put
  % insert only valid encoding vectors
  font_encoding_vector length 256 eq { newfont /Encoding font_encoding_vector put } if
  newfontname newfont definefont pop
} bind def

/Arial /font_Arial MF
%%EndProlog
%%Page: 1 1

/font_Arial findfont
10 scalefont
setfont
"""


# Pour la posterite: un debut de frontend Cairo

def cairoPrintHeader(landscape = False, format = "ps", filename = '/dev/stdout'):

	format = format.strip().lower()

	if landscape:
		width,height = 29.7, 21
	else:
		width,height = 21, 29.7

	if format == 'ps':
		surface = cairo.PSSurface(filename, width*72/2.54, height*72/2.54)
	elif format == 'pdf':
		surface = cairo.PSSurface(filename, width*72/2.54, height*72/2.54)
	else:
		print("Unknown output style '%s'" % format, file=sys.stderr)
		sys.exit(1)

	global ctx
	ctx = cairo.Context(surface)
	ctx.scale(width, height)
	ctx.set_line_width(0.001)

	initColor()


def cairoSetLineColor(C):
	global ctx

	if type(C) != tuple:
		C = getColor(C, (0,0,0))
	ctx.set_source_rgb(*C)

def cairoDrawLine(X, Y, L, H, C):

	global ctx

	setLineColor(C)
	ctx.move_to(X, Y)
	ctx.rel_line_to(L, H)
	ctx.stroke()


def cairoDrawBox(X, Y, L, H, Cb, Cr):

	setLineColor(Cb)
	ctx.rectangle(X, Y, L, H)

	if Cr != None:
		ctx.set_source(cairo.SolidPattern(Cr[0], Cr[1], Cr[2], 1))
		ctx.fill()

