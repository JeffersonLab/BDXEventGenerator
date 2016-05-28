#!/usr/local/bin/python

import string,math,os,sys,subprocess,shutil
import argparse

from CardsUtils import *
from LHEUtils import *
from DumpUtils import *

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

print "BDX EVENT GENERATOR PYTHON SCRIPT"

#check current working dir
EventGeneratorLocation  = os.environ['BDX_EVENT_GENERATOR']
MadGraphLocation        = EventGeneratorLocation+"/AprimeAlAlpha1"
MadGraphEventsLocation  = MadGraphLocation+"/Events"
MadGraphCardsLocation   = MadGraphLocation+"/Cards"
DetectorInteractionLocation   =  EventGeneratorLocation+"/DetectorInteraction"
DetectorInteractionEventsLocation = DetectorInteractionLocation+"/Events"

os.chdir(EventGeneratorLocation)


#no_showering=False

parser = argparse.ArgumentParser(description='BDX event generator')

parser.add_argument('--run_name',type=str,required=True,help="Run name",default='BDX');

parser.add_argument('--run_card',type=str, required=True,help='Run card to use')
parser.add_argument('--param_card',type=str,help='Param card to use',default='Cards/param_card.dat');
parser.add_argument('--proc_card',type=str,help='Proc card to use',default='Cards/proc_card.dat');
parser.add_argument('--max_attempts',type=int,default=3,help="Number of attempts per bin");
parser.set_defaults(no_showering=False);

parser.add_argument('--force_no_showering',dest='no_showering',action='store_true',help='Ignore showering effects, no matter what is used in the run_card');
args = parser.parse_args()

no_showering = args.no_showering;

run_name = args.run_name;
run_card_name=args.run_card;
param_card_name=args.param_card;
proc_card_name=args.proc_card;
max_attempts=args.max_attempts;

Energy = [];               #Array  of bin-center energy (GeV)
nGeneratedEvents = [];     #Array  of number of generated events per energy bin
nRequestedEvents = [];     #Array  of number of events that needs to be taken from this bin
Sigmas = [];               #Array  of total chi-chi production cross-section sigma(Ei) in pbarn per energy bin
Density= [];               #Array  of the quantity <dn/dE>*deltaE computed at the bin center, where <dn/dE> is the t-integrated, energy-differential, distribution of electrons in the dump per incident electron.
Weights= [];               #Array  of the quantity sigma(Ei)*<dn/dE>*deltaE computed at the bin center
NormWeights= [];           #Array  of the normalized weights (to the weights sum)
totalWeight=0;
NrequestedTOT=0;
NTOT = 0;
#Check if the cards exist
if (os.path.isfile(run_card_name)==False):
     print >> sys.stderr, bcolors.FAIL,"error: run card "+run_card_name+" does not exists!",bcolors.ENDC
     sys.exit(0)
     
if (os.path.isfile(proc_card_name)==False):
     print >> sys.stderr, bcolors.FAIL,"error: proc card "+proc_card_name+" does not exists!",bcolors.ENDC
     sys.exit(0)

if (os.path.isfile(param_card_name)==False):
     print >> sys.stderr,  bcolors.FAIL,"error: run_card "+param_card_name+" does not exists!",bcolors.ENDC
     sys.exit(0)
     
#Get the number of requested events from the run card
NrequestedTOT =  GetRequestedEvents(run_card_name)   
print bcolors.OKGREEN,"Requested: ",NrequestedTOT," events "

#next, check if, in the run_card, an option to use electron showering in the dump was implemented
#This will return a bool (yes/no), the Emin value,the Emax==Ebin value, the number of bins to use
UseElectronShowering,Emin,Emax,nbins = CheckElectronShowering(run_card_name,True)
deltaE = (Emax-Emin)/(nbins);

if (no_showering==True):
    print bcolors.WARNING,"Overriding electron showering, not using it",bcolors.ENDC
    UseElectronShowering=False


 
#Now, two opposite cases can happen. 
#First case: no showering case.
#This is the simplest scenario. We just need to call MadGraph once, with the RunCard provided
#1) Copy all the cards to the madgraph cards folder
#2) Run MadGraph
#3) Remove the Cards from the madgraph cards folder
if (UseElectronShowering==False):
    shutil.copy(run_card_name,MadGraphCardsLocation);
    shutil.copy(proc_card_name,MadGraphCardsLocation);
    shutil.copy(param_card_name,MadGraphCardsLocation);
    os.chdir(MadGraphLocation)
    print bcolors.OKGREEN,"Calling MadGraph once for run name: ",run_name,bcolors.ENDC
    command = "./bin/generate_events_new 0 "+run_name;
    os.system(command)
    command = "cd Events ; gzip -d "+run_name+"_unweighted_events.lhe.gz"     #Ok, this seems stupid since I am using gzip -d after MadGraph did a gzip. But I do not want to touch madgraph
    lhefname = MadGraphEventsLocation+"/"+run_name+"_unweighted_events.lhe"
    os.system(command)
    os.chdir(EventGeneratorLocation)
    
    os.remove(MadGraphCardsLocation+"/run_card.dat");
    os.remove(MadGraphCardsLocation+"/proc_card.dat");
    os.remove(MadGraphCardsLocation+"/param_card.dat");
    
    nGeneratedEventsThisRun,sigmaThisRun = GetGeneratedEventsNandSigma(lhefname)
    print bcolors.OKGREEN,"DONE: number of generated events is: ",nGeneratedEventsThisRun," cross section pb is: ",sigmaThisRun,bcolors.ENDC
    nGeneratedEvents.append(nGeneratedEventsThisRun);
    Sigmas.append(sigmaThisRun);
    Energy.append(Emax);
    Density.append(1.);  #It is correct to keep this to one -> There is exactly one electron per incident electron.
    Weights.append(sigmaThisRun);
    NormWeights.append(1.);
    print bcolors.OKGREEN,"LHE file was written",bcolors.ENDC  
#Second: showering case    
#In this case, one needs to call MadGraph more than once, by changing the primary beam energy 
else:
    print bcolors.OKGREEN,"Calling MadGraph ",nbins," times for run name: ",run_name,bcolors.ENDC
    
    for ii in range(nbins):

        flagThisLoopIteration=True;
        attemptsThisBin = 0;
        while(flagThisLoopIteration):
            Ei = Emin + deltaE/2 + deltaE * ii;                                                      #This is the energy to be used in this run
            CreateRunCardDifferentEnergy(Ei,run_card_name,MadGraphCardsLocation+"/run_card.dat");    #This function create the run card
            shutil.copy(proc_card_name,MadGraphCardsLocation);
            shutil.copy(param_card_name,MadGraphCardsLocation);
            print  bcolors.OKGREEN,"Energy: ",Ei,bcolors.ENDC
            os.chdir(MadGraphLocation)
            this_run_name=run_name+"_"+str(ii)
            command = "./bin/generate_events_new 0 "+this_run_name;  
            os.system(command)  
            command = "cd Events ; gzip -d "+this_run_name+"_unweighted_events.lhe.gz"     #Ok, this seems stupid since I am using gzip -d after MadGraph did a gzip. But I do not want to touch madgraph  
            lhefname = MadGraphEventsLocation+"/"+this_run_name+"_unweighted_events.lhe"
            os.system(command)
            
            os.chdir(EventGeneratorLocation)  
            os.remove(MadGraphCardsLocation+"/run_card.dat");
            os.remove(MadGraphCardsLocation+"/proc_card.dat");
            os.remove(MadGraphCardsLocation+"/param_card.dat");
            
            nGeneratedEventsThisRun,sigmaThisRun = GetGeneratedEventsNandSigma(lhefname)
            if (nGeneratedEventsThisRun==0):
                print bcolors.WARNING,"Error! 0 events generated. Trying again this bin. Next attempt is: ",str(attemptsThisBin),bcolors.ENDC
                attemptsThisBin = attemptsThisBin+1
                if (attemptsThisBin == max_attempts):
                    print bcolors.FAIL,"Error Error Error",bcolor.ENDC  
                    print bcolors.FAIL,"Error Error Error: max number of attempts exceeded. End Program here!",bcolors.ENDC
                    sys.exit(0)  
            else:   
                flagThisLoopIteration=False;  #This will break the while loop
                nGeneratedEvents.append(nGeneratedEventsThisRun);
                Energy.append(Ei);
                Sigmas.append(sigmaThisRun);
                DensityThisRun = dNdEIntegral(Ei)*deltaE;
                WeightThisRun = DensityThisRun * sigmaThisRun;
                Density.append(DensityThisRun);
                Weights.append(WeightThisRun);
                print bcolors.OKGREEN,"DONE. This attempt was: ",str(attemptsThisBin),bcolors.ENDC
                print bcolors.OKGREEN,"DONE: number of generated events is: ",nGeneratedEventsThisRun," cross section pb is: ",sigmaThisRun," density is: ",DensityThisRun,bcolors.ENDC
                print ""
        
    #At this point, we have generated madgraph events for all the energies Ei. We also have ALL the cross-sections and all the dn/dE|E=Ei    
 
    #Compute total weight: this is the total cross-section in pbarn per incident electron
    for ii in range(nbins):
        totalWeight +=Weights[ii];
   
    #Compute normalized weight per bin
    for ii in range(nbins):
        NormWeights.append(Weights[ii]/totalWeight)
    
    #Compute requested events per bin
    for ii in range(nbins):
         nRequestedEvents.append(int(NormWeights[ii]*NrequestedTOT));
         print bcolors.OKGREEN,"Bin ",ii," with energy ",Energy[ii],"sigma: ",Sigmas[ii]," density: ",Density[ii]," requested: ",nRequestedEvents[ii]," events, available: ",nGeneratedEvents[ii]," events ",bcolors.ENDC
     
    #Now check if ALL the bins have the necessary events. If not, we need to rescale.
    frac=1.
    fracMin=1.
    for ii in range(nbins):
        if (nRequestedEvents[ii]> nGeneratedEvents[ii]):
            frac=(1.*nGeneratedEvents[ii])/nRequestedEvents[ii];
            if (frac<fracMin):
                fracMin=frac;
    if (fracMin<1):
        print bcolors.FAIL,"Can generated only ",fracMin," of the requested events",bcolors.ENDC
        for ii in range(nbins):
            nRequestedEvents[ii]=int(nRequestedEvents[ii]*fracMin-1);  #The -1 is a conservative rounding

    for ii in range(nbins):
        NTOT+=nRequestedEvents[ii];
    print bcolors.OKGREEN,"TOTAL events written to LHE: ",NTOT,bcolors.ENDC   
    print bcolors.OKBLUE,"Mediated cross-section per incident electron: ",totalWeight," (pbarn) ",bcolors.ENDC     
    #Now, create the FINAL events file in LHE format, by keeping only the first nRequestedEvents from each file
    #First, create the LHE file, copying the header from the FIRST produced file
    lhefname = MadGraphEventsLocation+"/"+run_name+"_unweighted_events.lhe"
    lhefnameThisRun = MadGraphEventsLocation+"/"+run_name+"_0_unweighted_events.lhe"
    CreateLHEFileHeaderOnly(Emax,NTOT,totalWeight,lhefnameThisRun, lhefname)
    for ii in range(nbins):    
         lhefnameThisRun = MadGraphEventsLocation+"/"+run_name+"_"+str(ii)+"_unweighted_events.lhe"
         AppendEventsToLHEFileNewWeight(nRequestedEvents[ii],totalWeight/NTOT,lhefnameThisRun,lhefname)
    CloseLHEFile(lhefname)
    print bcolors.OKGREEN,"LHE file was written...",bcolors.ENDC         

#   At this point, we have generated chi-chi events in MadGraph, and we have also computed correctly the total cross-section.
    #Second part: detector interaction  
print bcolors.OKGREEN," Compiling DetectorInteraction ",bcolors.ENDC
os.chdir(DetectorInteractionLocation);
#Force recompilation
command = "make clean ; make"
os.system(command)
#run it 
print bcolors.OKGREEN," Running DetectorInteraction ",bcolors.ENDC
lhefname = MadGraphEventsLocation+"/"+run_name+"_unweighted_events.lhe"
lhefnameOUT = DetectorInteractionEventsLocation+"/"+run_name
command = "./DetectorInteraction.exe "+lhefname+" "+lhefnameOUT  
os.system(command)
os.chdir(EventGeneratorLocation)
