#!/usr/local/bin/python

import string,math,os,sys,subprocess,shutil
import os.path

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False




def createRunCard(nevents,eBeam,ldet,fiduciallx,fiducially,fiduciallz,procid):
        print "CREATE RUN CARD"
	runcardIN = open("run/run_card_base.dat", "r")
	runcardOUT = open("run/run_card.dat","w")
	lines = runcardIN.readlines()
	iii = 0
	for line in lines:
		x = line.split()
		if (len(x)>0):
			if (x[0][0]=="#"):
				runcardOUT.write(line)
			elif	(is_number(x[0])==False):
				runcardOUT.write(line)
			elif ((is_number(x[0])==True) and (x[1]=="=")):
		 		parName=x[2]	
				if (parName=="nevents"):
					x[0]=str(nevents)
				elif (parName=="ebeam1"):
					x[0]=str(eBeam)
				elif (parName=="ldet"):
					x[0]=str(ldet)
				elif (parName=="fiduciallx"):
					x[0]=str(fiduciallx)
				elif (parName=="fiducially"):
					x[0]=str(fiducially)
				elif (parName=="fiduciallz"):
					x[0]=str(fiduciallz)
				elif (parName=="procid"):
					x[0]=str(procid)	
				runcardOUT.write(" ")
				for word in x:
					runcardOUT.write(word+" ")
				runcardOUT.write("\n")
                        else:
                            paramcardOUT.write(line)
                        
	runcardIN.close()
	runcardOUT.close()


def createParamCard(aMass,chiMass,epsilon,alphaD):
        print "CREATE PARAM CARD"
	paramcardIN = open("run/param_card_base.dat", "r")
	paramcardOUT = open("run/param_card.dat","w")
	lines = paramcardIN.readlines()
	iii = 0
	for line in lines:
		x = line.split()
		if (len(x)>0):
			if (x[0][0]=="#"):
				paramcardOUT.write(line)
			elif	(is_number(x[0])==False):
				paramcardOUT.write(line)
			elif ((is_number(x[0])==True) and (is_number(x[1])==True) and (x[2]=="#")):
				parName=x[3]
				if (parName=="FMASS"):
					x[1]=str(chiMass)
				elif (parName=="APMASS"):
					x[1]=str(aMass)
				elif (parName=="epsilon"):
					x[1]=str(epsilon)
				elif (parName=="alphaD"):
					x[1]=str(alphaD)
				paramcardOUT.write(" ")
				for word in x:
					paramcardOUT.write(word+"  ")
				paramcardOUT.write("\n")
                        else:
                            paramcardOUT.write(line)
	paramcardIN.close()
	paramcardOUT.close()

def createCards(fname):
    print "CREATE CARDS"
    print fname
    if (os.path.isfile(fname)==False):
        sys.exit(1)        
    else:
        cardIN = open(fname, "r")
       	lines = cardIN.readlines() 
        if (len(lines)!=1):
            sys.exit(2)
        line=lines[0]
        x = line.split()
        createRunCard(x[5],x[6],x[7],x[8],x[9],x[10],x[11])
        createParamCard(x[1],x[2],x[3],x[4])
    sys.exit(0)
