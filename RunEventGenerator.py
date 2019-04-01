#!/usr/local/bin/python

import string,math,os,sys,subprocess,shutil
import argparse

from CardsUtils import *
import LHEUtils 
#from DumpUtils import *
from PositronAnnihilationUtils import *
from ROOT import *



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


parser = argparse.ArgumentParser(description='BDX event generator')

parser.add_argument('--run_name',type=str,required=True,help="Run name",default='BDX');
parser.add_argument('--run_card',type=str, required=True,help='Run card to use')
parser.add_argument('--param_card',type=str,help='Param card to use',default='Cards/param_card.dat');
parser.add_argument('--proc_card',type=str,help='Proc card to use',default='Cards/proc_card.dat');
parser.add_argument('--max_attempts',type=int,default=3,help="Number of attempts per bin");
parser.set_defaults(no_showering=False);
parser.add_argument('--force_no_showering',dest='no_showering',action='store_true',help='Ignore showering effects, no matter what is used in the run_card');
parser.add_argument('--root_file',type=str,required=True,help="The root file with the data from the beam-dump simulation. In particular, the histogram (named hEall) containing dN/dE per incident electron and the histogram (named hE_angle_all) with the dN/dEdcosTheta. These two MUST have same binning wrt to x axis (energy axis)")
parser.add_argument('--add_annihilation',dest="annihilation",action='store_true',help='Add also the contribution from e+e- -> Aprime -> chi chiBar if set. It requires to run the showering code')
args = parser.parse_args()

no_showering = args.no_showering;
annihilation = args.annihilation;

run_name = args.run_name;
run_card_name=args.run_card;
param_card_name=args.param_card;
proc_card_name=args.proc_card;
max_attempts=args.max_attempts;

#Write in a global var the Mchi,Ma,eps,alphaD
mN = 0; #Target nucleous mass
mChi=0;
mA=0;
eps=0;
alphaD=0;
 
Energy = [];               #Array  of bin-center energy (GeV)
nGeneratedEvents = [];     #Array  of number of generated events per energy bin
nIdeallyTotalEvents = [];  #Array  of number of total events corresponding to this bin (nGeneratedEvents / normWeight)
nRequestedEvents = [];     #Array  of number of events that needs to be taken from this bin
Sigmas = [];               #Array  of total chi-chi production cross-section sigma(Ei) in pbarn per energy bin
Density= [];               #Array  of the quantity <dn/dE>* computed at the bin center, where <dn/dE> is the t-integrated, energy-differential, distribution of electrons in the dump per incident electron.
Weights= [];               #Array  of the quantity sigma(Ei)*<dn/dE>*deltaE computed at the bin center
NormWeights= [];           #Array  of the normalized weights (to the weights sum)
totalWeight=0;
NrequestedTOT=0;
NrequestedMAX=0;
NTOT = 0;

#A.C. 25/3/2019: quantities for the annihilation
WeightsAnnihilation=[];
totalWeightAnnihilation=0;
nRequestedEventsAnnihilation = [];
NTOTAnnihilation = 0;

deltaE=0;

#root stuff
rootFile=0;
hEall=0;
hE_angle_all=0;
hEallP=0;
hE_angle_allP=0;

BigNumber=9999999999;

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
print bcolors.OKGREEN,"Requested: ",NrequestedTOT," events ",bcolors.ENDC

#next, check if, in the run_card, an option to use electron showering in the dump was implemented
#This will return a bool (yes/no) and the BeamEnergy
UseElectronShowering,Ebeam = CheckElectronShoweringNew(run_card_name,True)


if (no_showering==True):
    print bcolors.WARNING,"Overriding electron showering, not using it",bcolors.ENDC
    UseElectronShowering=False
    if (annihilation==True):
        print bcolors.FAIL,"ERROR, can't compute annihilation when not using showering. End",bcolors.ENDC
        exit
#Get from the provided root files the histogram with the <dN/dE> - integrated over t - , normalized to incident electrons
else:
    rootFile=TFile(args.root_file);
    hEall=rootFile.Get("hEall");
    hEallP=rootFile.Get("hEallP");
    hE_angle_all=rootFile.Get("hE_angle_all")
    hE_angle_allP=rootFile.Get("hE_angle_allP")
    deltaE = hEall.GetXaxis().GetBinWidth(1) #this is the bin width
    nbins = hEall.GetNbinsX(); #this is the number of bins
    Emin= hEall.GetBinCenter(1)-deltaE/2;
    Emax= hEall.GetBinCenter(nbins)+deltaE/2;
    
    #Check Emax
    if (Ebeam!=Emax):
        print bcolors.WARNING," Ebeam in the run card: ",Ebeam," different from Emax in the histogram: ",Emax,bcolors.ENDC
        print bcolors.WARNING," using as Emax the one from histogram ",bcolors.ENDC
    #need to compute here the number of bins
    print "deltaE is: ",deltaE," nBins is: ",nbins, " eMin is: ",Emin," eMax is: ",Emax


#Now, two opposite cases can happen. 
#First case: no showering case.
#This is the simplest scenario. We just need to call MadGraph once, with the RunCard provided
#1) Copy all the cards to the madgraph cards folder
#2) Run MadGraph
#3) Remove the Cards from the madgraph cards folder
if (UseElectronShowering==False):
    shutil.copy(run_card_name,MadGraphCardsLocation+"/run_card.dat");
    shutil.copy(proc_card_name,MadGraphCardsLocation+"/proc_card.dat");
    shutil.copy(param_card_name,MadGraphCardsLocation+"/param_card.dat");
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
    
    nGeneratedEventsThisRun,sigmaThisRun = LHEUtils.GetGeneratedEventsNandSigma(lhefname)
    print bcolors.OKGREEN,"DONE: number of generated events is: ",nGeneratedEventsThisRun," cross section pb is: ",sigmaThisRun,bcolors.ENDC
    nGeneratedEvents.append(nGeneratedEventsThisRun);
    Sigmas.append(sigmaThisRun);
    Energy.append(Ebeam);
    Density.append(1.);  #It is correct to keep this to one -> There is exactly one electron per incident electron.
    Weights.append(sigmaThisRun);
    NormWeights.append(1.);
    print bcolors.OKGREEN,"LHE file was written",bcolors.ENDC  
#Second: showering case    
#In this case, one needs to call MadGraph more than once, by changing the primary beam energy 
else:
    print bcolors.OKGREEN,"Calling MadGraph ",nbins," times for run name: ",run_name,bcolors.ENDC
    
    #check the kinematics
    mChi = GetChiMass(param_card_name);
    mA   = GetAprimeMass(param_card_name);
    mN   = GetTargetMass(run_card_name);
    alphaD = GetAlphaDark(param_card_name);
    eps    = GetEpsilon(param_card_name);
    Znuc   = GetZnuc(param_card_name);
    print bcolors.OKGREEN,"mCHI: ",mChi,"mA: ",mA,"mN: ",mN,"alphaD: ",alphaD,"eps: ","Znuc: ",Znuc,eps,bcolors.ENDC
    
    
    for ii in range(nbins):

        flagThisLoopIteration=True;
        attemptsThisBin = 0;
        while(flagThisLoopIteration):
#            Ei = Emin + deltaE/2 + deltaE * ii;        #This is the energy to be used in this run
            Ei = hEall.GetBinCenter(ii+1);
            print  bcolors.OKGREEN,"Energy: ",Ei,bcolors.ENDC
            CreateRunCardDifferentEnergy(Ei,run_card_name,MadGraphCardsLocation+"/run_card.dat");    #This function create the run card
            shutil.copy(proc_card_name,MadGraphCardsLocation);
            shutil.copy(param_card_name,MadGraphCardsLocation);
            os.chdir(MadGraphLocation)
            this_run_name=run_name+"_"+str(ii)
            Ethr = 2*mChi*(1+mChi/mN);
            print bcolors.OKGREEN,"Chi mass: ",mChi," N mass:",mN,"thr: ",Ethr,bcolors.ENDC
            if (Ethr > Ei):
                Energy.append(Ei);
                DensityThisRun = hEall.GetBinContent(hEall.FindBin(Ei))*deltaE;    #This is for histogram
                nGeneratedEvents.append(0); #No events this run
                Sigmas.append(0.);          #No events this run
                Density.append(DensityThisRun);
                WeightThisRun =0.;
                Weights.append(WeightThisRun);
                print bcolors.WARNING,"Kinematic constraint for Ei: ",Ei," not satisfied, emin: ",Ethr,bcolors.ENDC
                os.chdir(EventGeneratorLocation)  
                os.remove(MadGraphCardsLocation+"/run_card.dat");
                os.remove(MadGraphCardsLocation+"/proc_card.dat");
                os.remove(MadGraphCardsLocation+"/param_card.dat");
                flagThisLoopIteration=False;  #This will break the while loop
                continue;#This will break the while loop

            command = "./bin/generate_events_new 0 "+this_run_name;  
            os.system(command)  
            command = "cd Events ; gzip -d "+this_run_name+"_unweighted_events.lhe.gz"     #Ok, this seems stupid since I am using gzip -d after MadGraph did a gzip. But I do not want to touch madgraph  
            lhefname = MadGraphEventsLocation+"/"+this_run_name+"_unweighted_events.lhe"
            os.system(command)
            
            os.chdir(EventGeneratorLocation)  
            os.remove(MadGraphCardsLocation+"/run_card.dat");
            os.remove(MadGraphCardsLocation+"/proc_card.dat");
            os.remove(MadGraphCardsLocation+"/param_card.dat");     
            nGeneratedEventsThisRun,sigmaThisRun = LHEUtils.GetGeneratedEventsNandSigma(lhefname);
            if (nGeneratedEventsThisRun==0):
                print bcolors.WARNING,"Error! 0 events generated. Trying again this bin. Next attempt is: ",str(attemptsThisBin),bcolors.ENDC
                attemptsThisBin = attemptsThisBin+1
                if (attemptsThisBin == max_attempts):
                    print bcolors.FAIL,"Error Error Error",bcolors.ENDC  
                    print bcolors.FAIL,"Error Error Error: max number of attempts exceeded. End Program here!",bcolors.ENDC
                    sys.exit(0)  
            else:   
                flagThisLoopIteration=False;  #This will break the while loop
                nGeneratedEvents.append(nGeneratedEventsThisRun);
                print "This run generated ",nGeneratedEventsThisRun
                Energy.append(Ei);
                Sigmas.append(sigmaThisRun);
                DensityThisRun =  hEall.GetBinContent(hEall.FindBin(Ei))*deltaE;    #This is for histogram
          ##      DensityThisRun = dNdEIntegral(Ei)*deltaE;                         #This is for analytical
                WeightThisRun = DensityThisRun * sigmaThisRun;
                Density.append(DensityThisRun);
                Weights.append(WeightThisRun);

                #Rotate the events
                print bcolors.OKGREEN,"Now rotate events. Projecting bin: ",ii," Check, bin center is: ",hE_angle_all.GetXaxis().GetBinCenter(ii+1),bcolors.ENDC
                hAngleTMP=hE_angle_all.ProjectionY("_py",ii+1,ii+1)
                print lhefname
                LHEUtils.RotateLHEFEvents(lhefname,hAngleTMP)

                print bcolors.OKGREEN,"DONE. This attempt was: ",str(attemptsThisBin),bcolors.ENDC
                print bcolors.OKGREEN,"DONE: number of generated events is: ",nGeneratedEventsThisRun," cross section pb is: ",sigmaThisRun," density is: ",DensityThisRun,bcolors.ENDC
                print ""
        
    #At this point, we have generated madgraph events for all the energies Ei. We also have ALL the cross-sections and all the dn/dE|E=Ei    
    
    #A.C. on 25/3/2019: adding e+ e- -> A' -> chi chiBar
    if (annihilation==True):
        annihilationHandler=PositronAnnihilationSpectra(mA,mChi,eps,alphaD,hEallP)
        for ii in range(nbins):
            WeightsAnnihilation.append(annihilationHandler.xsectionAvg[ii]*Znuc); #this is the annihilation xsection integrated over T(E) for this bin. IMPORTANT: this is multiplied by Znuc - since the code handles normalization has Nav/A for the Breem. coherent process.
            totalWeightAnnihilation+=WeightsAnnihilation[ii];  #This is the annihilation xsection integrated over T(E) for all bins
        print bcolors.OKGREEN,"TOTAL Int(sigmaE*T(E))= ",totalWeightAnnihilation," (pbarn)",bcolors.ENDC
            
            
    #Compute total weight: this is the total cross-section in pbarn per incident electron
    for ii in range(nbins):
        totalWeight +=Weights[ii];

    #Compute normalized weight per bin
    for ii in range(nbins):
        NormWeights.append(Weights[ii]/totalWeight)
    
    #Compute requested events per bin. Actually, I want to maximize the number of events!!
    for ii in range(nbins):
        if (Weights[ii]<=0.): #case of bins with no xsection due to threshold
            nIdeallyTotalEvents.append(BigNumber)
            print bcolors.OKBLUE,"Bin ",ii," nIdeallyTotalevents: ",BigNumber," BIG NUMBER TO AVOID bin",bcolors.ENDC            
        else:
            nIdeallyTotalEvents.append(nGeneratedEvents[ii]/NormWeights[ii])
            print bcolors.OKBLUE,"Bin ",ii," nIdeallyTotalevents: ",nIdeallyTotalEvents[ii],bcolors.ENDC            
    nRequestedMAX=min(nIdeallyTotalEvents)
    nRequestedMAX=(int)(nRequestedMAX-1) #just for safety
    print bcolors.OKGREEN,"Min is: ",nRequestedMAX
    
    for ii in range (nbins):   
         nRequestedEvents.append(int(NormWeights[ii]*nRequestedMAX));
         print bcolors.OKGREEN,"Bin ",ii," with energy ",Energy[ii],"sigma: ",Sigmas[ii],"normweight: ",NormWeights[ii]," density: ",Density[ii]," requested: ",nRequestedEvents[ii]," events, available: ",nGeneratedEvents[ii]," events ",bcolors.ENDC
     
    #Now check if ALL the bins have the necessary events (should always be!). If not, we need to rescale.
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

    #At this point, I know the number of requested events per each bin.
    for ii in range(nbins):
        NTOT+=nRequestedEvents[ii];
    print bcolors.OKGREEN,"TOTAL Bremmstrahlung events written to LHE: ",NTOT,bcolors.ENDC   
    print bcolors.OKBLUE,"Mediated Bremmstrahlung cross-section per incident electron: ",totalWeight," (pbarn) ",bcolors.ENDC     
    
     #A.C. 25/3/2019: compute the number of requested events per bin for the annihilation
    if (annihilation == True):
        for ii in range(nbins):
            if (Weights[ii]>0):
                nRequestedEventsAnnihilation.append(int(nRequestedEvents[ii]*WeightsAnnihilation[ii]/Weights[ii]))
            else:
                nRequestedEventsAnnihilation.append(0)
            NTOTAnnihilation+=nRequestedEventsAnnihilation[ii];
            print bcolors.OKGREEN,"Bin ",ii," annihilation required: ",nRequestedEventsAnnihilation[ii]
            
        print bcolors.OKGREEN,"TOTAL Annihilation events written to LHE: ",NTOTAnnihilation,bcolors.ENDC   
        print bcolors.OKBLUE,"Mediated Annihilation cross-section per incident electron: ",totalWeightAnnihilation," (pbarn) ",bcolors.ENDC     
      
        #Very important: if annihilation is enabled, sum the corresponding total weight and the number of events
        NTOT = NTOT + NTOTAnnihilation;
        totalWeight = totalWeight + totalWeightAnnihilation
    
    print bcolors.OKGREEN,"TOTAL events written to LHE: ",NTOT,bcolors.ENDC   
    print bcolors.OKBLUE,"Mediated cross-section per incident electron: ",totalWeight," (pbarn) ",bcolors.ENDC     
    
    #Now, create the FINAL events file in LHE format, by keeping only the first nRequestedEvents from each file
    #First, create the LHE file, copying the header from the FIRST produced file
    lhefname = MadGraphEventsLocation+"/"+run_name+"_unweighted_events.lhe"
    lhefnameThisRun="";
    for ii in range(nbins):    
        if (nRequestedEvents[ii]>0):
            lhefnameThisRun = MadGraphEventsLocation+"/"+run_name+"_"+str(ii)+"_unweighted_events.lhe"
            break;
    if (lhefnameThisRun==""):
        print bcolors.FAIL,"Error can't get a good header! EXIT",bcolors.ENDC
        exit(1);
    LHEUtils.CreateLHEFileHeaderOnly(Emax,NTOT,totalWeight,lhefnameThisRun, lhefname)
    for ii in range(nbins):    
        if (nRequestedEvents[ii]>0):
            lhefnameThisRun = MadGraphEventsLocation+"/"+run_name+"_"+str(ii)+"_unweighted_events.lhe"
            LHEUtils.AppendEventsToLHEFileNewWeight(nRequestedEvents[ii],totalWeight/NTOT,lhefnameThisRun,lhefname)
        if (annihilation == True):
            if (nRequestedEventsAnnihilation[ii]>0):   
                hAngleTMP=hE_angle_allP.ProjectionY("_py",ii+1,ii+1)
                annihilationHandler.generateAndWriteEventsWRotation(nRequestedEventsAnnihilation[ii],ii,totalWeight/NTOT,hAngleTMP,lhefname)
            
    LHEUtils.CloseLHEFile(lhefname)
    print bcolors.OKGREEN,"LHE file was written...",bcolors.ENDC         

#    At this point, we have generated chi-chi events in MadGraph, and we have also computed correctly the total cross-section.
#    Second part: detector interaction  
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
