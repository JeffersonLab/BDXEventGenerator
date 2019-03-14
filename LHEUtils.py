#!/usr/local/bin/python
import string,math,os,sys
from ROOT import *

#This function reads the LHE file given and returns:
# the number of actually generated events present in the file
# the TOTAL cross-section for these events in pbarn 
def GetGeneratedEventsNandSigma(fname):
    nEvents = 0
    isMGGenerationInfo = False
    lhefile=open(fname,"r")
    
    lines =  lhefile.readlines()
    
    for line in lines:
        if (line.find("<MGGenerationInfo>") != -1):
            isMGGenerationInfo = True;
            continue;
        elif (line.find("</MGGenerationInfo>") != -1):
            isMGGenerationInfo = False;
            break;
        elif ((line.find(" Number of Events ") != -1)and(isMGGenerationInfo==True)):   #Line is "# Number of Events    :   N"
            x=line.split();
            nEvents=int(x[5]);
            continue;
        elif ((line.find(" Integrated weight (pb) ") != -1)and(isMGGenerationInfo==True)):   #Line is "#  Integrated weight (pb)  :  sigma"
            x=line.split();
            sigma=float(x[5]);
    return nEvents,sigma;

#This function creates a LHE file (dest_lhe_filename), copying just the header from the src_lhe_filename,
#but making the following substitutions:
##1) ebeam1 is changed to ebeam in the banner containing the run card
##2) In the section <MGGenerationInfo>, the number of events ant the integrated weight are changed to Nevents and sigma(in pbarn).
##3) In the section <MGGenerationInfo>, the truncated weight is put to 0
##4) In the section <MGGenerationInfo>, the unit wg is corrected to sigma/Nevents

def CreateLHEFileHeaderOnly(E0,Nevents,sigma,src_lhe_filename,dest_lhe_filename):
    src_lhefile=open(src_lhe_filename,"r");
    dest_lhefile=open(dest_lhe_filename,"w");
    isMGGenerationInfo = False
    lines =  src_lhefile.readlines()
    for line in lines:
        if (line.find(" ebeam1 ")!= -1):
            newline = " "+str(E0)+" = ebeam1 ! beam 1 energy in GeV\n"
            dest_lhefile.write(newline)
            continue
        elif (line.find("<MGGenerationInfo>") != -1):
            isMGGenerationInfo = True;
            dest_lhefile.write(line);    
            continue;
        elif (line.find("</MGGenerationInfo>") != -1):    
            isMGGenerationInfo = False;
            dest_lhefile.write(line);
            continue;   
        elif (line.find("</init>") != -1): 
            dest_lhefile.write(line);        
            break;
        elif ((line.find(" Number of Events ") != -1)and(isMGGenerationInfo==True)):   #Line is "# Number of Events    :   N"
            newline = "#  Number of Events        :        "+str(Nevents)+"\n"
            dest_lhefile.write(newline)
            continue;
        elif ((line.find(" Integrated weight (pb) ") != -1)and(isMGGenerationInfo==True)):   #Line is "#  Integrated weight (pb)  :  sigma"
            newline = "#  Integrated weight (pb)  :  "+str(sigma)+"\n"
            dest_lhefile.write(newline)
            continue;
        elif ((line.find(" Truncated wgt (pb) ") != -1)and(isMGGenerationInfo==True)):
            newline = "#  Truncated wgt (pb)  :  0\n"
            dest_lhefile.write(newline)  
            continue;
        elif ((line.find(" Unit wgt ") != -1)and(isMGGenerationInfo==True)):
            newline = "#  Unit wgt                :  "+str(sigma/Nevents)+"\n"
            dest_lhefile.write(newline);
            continue
        else:
            dest_lhefile.write(line);
            continue;
    src_lhefile.close();
    dest_lhefile.close();
 
#The following function appends to dest_lhe_filename the first Nevents of src_lhe_filename.
#It is basically a copy-to function    
def AppendEventsToLHEFile(Nevents,src_lhe_filename,dest_lhe_filename):    
    src_lhefile=open(src_lhe_filename,"r");
    dest_lhefile=open(dest_lhe_filename,"a");
    isEvent = False
    NwrittenEvents=0;
    
    lines =  src_lhefile.readlines()
    
    
    for line in lines:
        if (line.find("<event>") != -1):
            isEvent = True
        elif (line.find("</event>") != -1):
            isEvent = False
            dest_lhefile.write(line)
            NwrittenEvents=NwrittenEvents+1;
        if (isEvent==True):
            dest_lhefile.write(line)
        if (NwrittenEvents==Nevents):
            break;

#The following function appends to dest_lhe_filename the first Nevents of src_lhe_filename.
#For EACH event, the event weight is modified, and set to weight
def AppendEventsToLHEFileNewWeight(Nevents,weight,src_lhe_filename,dest_lhe_filename):    
    src_lhefile=open(src_lhe_filename,"r");
    dest_lhefile=open(dest_lhe_filename,"a");
    isEvent = False
    isHeader = False
    NwrittenEvents=0;
    lines =  src_lhefile.readlines()
    for line in lines:
        if (line.find("<event>") != -1):
            isEvent = True
            isHeader = True
        elif (isHeader):
            isHeader=False;
            x=line.split()              #This is the header line. Nparticles 0 eventWeight scale alphaEM alphaS
            newline=x[0]+" "+x[1]+" "+str(weight)+" "+x[3]+" "+x[4]+" "+x[5]+"\n";
            dest_lhefile.write(newline)
            continue #So we do not write this twice
        elif (line.find("</event>") != -1):
            isEvent = False
            dest_lhefile.write(line)
            NwrittenEvents=NwrittenEvents+1;
        if (isEvent==True):
            dest_lhefile.write(line)
        if (NwrittenEvents==Nevents):
            break;        

def RotateLHEFEvents(src_lhe_filename,hAngleTMP):
    src_lhefile=open(src_lhe_filename,"r");
    dest_lhefile=open(src_lhe_filename+".tmp","a");
    isFirstParticle = False
    isEvent = False
    isHeader = False
    lines =  src_lhefile.readlines()
    random=TRandom3(0)
    for line in lines:
        if (line.find("<event>") != -1):
            isFirstParticle = True
            isEvent = True
            isHeader = True
            dest_lhefile.write(line)
        elif (isHeader):
            isHeader=False;
            dest_lhefile.write(line)
        elif (line.find("</event>") != -1):    
            isEvent = False
            dest_lhefile.write(line)
        elif (isEvent):
            x=line.split()  
            px=(float)(x[6])
            py=(float)(x[7])
            pz=(float)(x[8])
            vector=TVector3(px,py,pz)
            if (isFirstParticle):
                ctheta=hAngleTMP.GetRandom()
                theta=TMath.ACos(ctheta)           #This returns 0 if ctheta > 1
                phi=random.Uniform(0,2*TMath.Pi())
                isFirstParticle=False
            vector.RotateY(theta)
            vector.RotateZ(phi)
            newline=x[0]+" "+x[1]+" "+x[2]+" "+x[3]+" "+x[4]+" "+x[5]+" "+str(vector.Px())+" "+str(vector.Py())+" "+str(vector.Pz())+" "+x[9]+" "+x[10]+" "+x[11]+" "+x[12]+"\n";
            dest_lhefile.write(newline)
            continue;
        else:
            dest_lhefile.write(line)
    src_lhefile.close()
    dest_lhefile.close()
    os.rename(src_lhe_filename+".tmp",src_lhe_filename)

def CloseLHEFile(dest_lhe_filename):        
    dest_lhefile=open(dest_lhe_filename,"a");
    line="</LesHouchesEvents>"
    dest_lhefile.write(line)
    dest_lhefile.close()         
