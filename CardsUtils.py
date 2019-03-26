#!/usr/local/bin/python
import string,math,os,sys


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def GetRequestedEvents(run_card_name):
    name = ""
    val = 0
    nevt = 1000;
    
    run_card=open(run_card_name,"r")
    lines = run_card.readlines()
    for line in lines:
        x = line.split()
        if (len(x)==0):
             continue;
        elif (is_number(x[0])and(len(x)>=2)): #This line reads a number
            val = x[0];
            name = x[2];
            if (name=="nevents"):
                nevt = int(val)
                break
    return nevt

def GetTargetMass(run_card_name,verbose=False):
    name = ""
    run_card=open(run_card_name,"r")
    lines = run_card.readlines()
    for line in lines:
        x = line.split()
        if (len(x)==0):
            continue;
        if (is_number(x[0])and(len(x)>=2)): #This line reads a number
            val = x[0];
            name = x[2];
        if (name=="mbeam2"):
            M = float(val);
            return M
    if (verbose):
        print "Error in GetTargetMass"
        return 0.

    
def GetChiMass(param_card_name,verbose=False):
    name = ""
    param_card=open(param_card_name,"r")
    lines = param_card.readlines();
    for line in lines:
        x = line.split()
        if (len(x)==0):
            continue;
        if ((is_number(x[0])==True) and (is_number(x[1])==True) and (x[2]=="#")):
            parName=x[3]
            if (parName=="FMASS"):
                M=float(x[1])
                return M
    if (verbose):
        print "Error in getChiMass"
    return 0.

def GetAprimeMass(param_card_name,verbose=False):
    name = ""
    param_card=open(param_card_name,"r")
    lines = param_card.readlines();
    for line in lines:
        x = line.split()
        if (len(x)==0):
            continue;
        if ((is_number(x[0])==True) and (is_number(x[1])==True) and (x[2]=="#")):
            parName=x[3]
            if (parName=="APMASS"):
                M=float(x[1])
                return M
    if (verbose):
        print "Error in getAprimeMass"
    return 0.

def GetAlphaDark(param_card_name,verbose=False):
    name = ""
    param_card=open(param_card_name,"r")
    lines = param_card.readlines();
    for line in lines:
        x = line.split()
        if (len(x)==0):
            continue;
        if ((is_number(x[0])==True) and (is_number(x[1])==True) and (x[2]=="#")):
            parName=x[3]
            if (parName=="alphaD"):
                alphaD=float(x[1])
                return alphaD
    if (verbose):
        print "Error in GetAlphaDark"
    return 0.

def GetEpsilon(param_card_name,verbose=False):
    name = ""
    param_card=open(param_card_name,"r")
    lines = param_card.readlines();
    for line in lines:
        x = line.split()
        if (len(x)==0):
            continue;
        if ((is_number(x[0])==True) and (is_number(x[1])==True) and (x[2]=="#")):
            parName=x[3]
            if (parName=="epsilon"):
                eps=float(x[1])
                return eps
    if (verbose):
        print "Error in GetEpsilon"
    return 0.

#This function reads the provided run_card, and checks if the electron showering option was asked for - or not.
#The relevant informations are in the "<BDX>" section
#This will return:
# useShowering: bool,   use or not the showering
# Emin        : double, minimum electron energy
# Emax        : double, maximum electron energy
# n           : number of bins in the (Emin,Emax) interval 
def CheckElectronShowering(run_card_name,verbose=False):
    
    isBDX = False;
  
    name = ""
    val = 0
    
    useShowering = False;
    Emax=11;
    Emin=1;
    n = 10;
    
    run_card=open(run_card_name,"r")
    lines = run_card.readlines()
    
    for line in lines:
        x = line.split()
        if (len(x)==0):
             continue;
        elif (line.find("<BDX>") != -1):  #This line starts the BDX part
            isBDX=True
        
        elif (line.find("</BDX>") != -1): #This line ends the BDX part
            isBDX=False

        elif (is_number(x[0])and(len(x)>=2)): #This line reads a number
            val = x[0];
            name = x[2];
            
        #Now check values
        #E0
        if (name=="ebeam1"):
            Emax = float(val);
        elif ((isBDX==True)and(name=="USESHOWERING")):
            if (int(val)>=1):
                 useShowering = True;
            else:
                 useShowering = False;
        elif ((isBDX==True)and(name=="EMINSHOWERING")):
            Emin = float(val);
        elif ((isBDX==True)and(name=="NBINSSHOWERING")):
            n = int(val);
    run_card.close();
    
    if (verbose==True):
        print "CheckElectronShowering:"
        print "useShowering: ",useShowering
        print "Ebeam:        ",Emax
    
    return useShowering,Emax,n
    

#This function reads the provided run_card, and checks if the electron showering option was asked for - or not.
#The relevant informations are in the "<BDX>" section
#This will return:
# useShowering: bool,   use or not the showering
# Ebeam        : double, beam energy
def CheckElectronShoweringNew(run_card_name,verbose=False):
    
    isBDX = False;
  
    name = ""
    val = 0
    
    useShowering = False;
    Ebeam=11;
  
    
    run_card=open(run_card_name,"r")
    lines = run_card.readlines()
    
    for line in lines:
        x = line.split()
        if (len(x)==0):
             continue;
        elif (line.find("<BDX>") != -1):  #This line starts the BDX part
            isBDX=True
        
        elif (line.find("</BDX>") != -1): #This line ends the BDX part
            isBDX=False

        elif (is_number(x[0])and(len(x)>=2)): #This line reads a number
            val = x[0];
            name = x[2];
            
        #Now check values
        #E0
        if (name=="ebeam1"):
            Ebeam = float(val);
        elif ((isBDX==True)and(name=="USESHOWERING")):
            if (int(val)>=1):
                 useShowering = True;
            else:
                 useShowering = False;
    run_card.close();
    
    if (verbose==True):
        print "CheckElectronShowering:"
        print "useShowering set to: ",useShowering
        print "Ebeam:        ",Ebeam
    
    return useShowering,Ebeam



def CreateRunCardDifferentEnergy(Ei,src_run_card_name,dest_run_card_name):    
    src_run_card=open(src_run_card_name,"r")
    dest_run_card=open(dest_run_card_name,"w")
    lines = src_run_card.readlines()
    name = ""
    for line in lines:
        x = line.split()
        if (len(x)==0):
             continue;
        elif ((is_number(x[0])and(len(x)>=2))): #This line reads a number
            val = x[0];
            name = x[2];
        if (name == "ebeam1"):
            newline = " "+str(Ei)+" = ebeam1 ! beam 1 energy in GeV\n"
            dest_run_card.write(newline)
        else:
            dest_run_card.write(line)
            
    src_run_card.close()
    dest_run_card.close()
    
    
def createParamCard(aMass,chiMass,epsilon,alphaD):
    print "CREATE PARAM CARD"
    paramcardIN = open("Cards/param_card.dat", "r")
    paramcardOUT = open("run/param_card.dat","w")
    lines = paramcardIN.readlines()
    iii = 0
    for line in lines:
        x = line.split()
        if (len(x)>0):
            if (x[0][0]=="#"):
                paramcardOUT.write(line)
            elif    (is_number(x[0])==False):
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
    
    
def createRunCard(nevents,eBeam,ldet,fiduciallx,fiducially,fiduciallz,procid):
    print "CREATE RUN CARD"
    runcardIN = open("Cards/run_card.dat", "r")
    runcardOUT = open("run/run_card.dat","w")
    lines = runcardIN.readlines()
    iii = 0
    for line in lines:
        x = line.split()
        if (len(x)>0):
            if (x[0][0]=="#"):
                runcardOUT.write(line)
            elif    (is_number(x[0])==False):
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
    
def createCards(fname):
    print "CREATE CARDS"
    print fname
    if (os.path.isfile(fname)==False):
        print "file does not exist"
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
    
    
