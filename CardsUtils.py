#!/usr/local/bin/python
import string,math,os
from utils import is_number



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
        print "Emin:         ",Emin
        print "Emax:         ",Emax
        print "n:            ",n
    
    return useShowering,Emin,Emax,n
    
    
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
    
    
    