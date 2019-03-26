#!/usr/local/bin/python
import string,math,os,sys
from ROOT import *

class mbcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

#All masses in GeV
#EA,EB is the lower-upper limit for positron energies in GeV in lab frame
class PositronAnnihilationKin:
    Me=0.511E-3;
    alpha=1./137;
    def __init__(self,Ma,Mchi,eps,alphaD,EA,EB):
        self.Ma=Ma
        self.Mchi=Mchi
        self.eps=eps
        self.alphaD=alphaD
        self.EA=EA
        self.EB=EB
        #Compute A' width to chi chiBar        
        if (1-4.*self.Mchi*self.Mchi/(self.Ma*self.Ma) < 0):
            self.widthI=0
        else:
            self.widthI=self.alphaD*self.Ma/3.*math.sqrt(1-4.*self.Mchi*self.Mchi/(self.Ma*self.Ma))*(1+2.*self.Mchi*self.Mchi/(self.Ma*self.Ma));

        self.widthV=self.alpha*self.eps*self.eps*self.Ma/3.;
        self.width=self.widthI+self.widthV;
        
        self.sdXsec=2*((EA+EB)/2.)*PositronAnnihilationKin.Me
        
        #Here the variable is the positron energy in lab frame in GeV
        self.xsecF = TF1("PositronAnnihilationTotXsection",self.getXsectionE,EA,EB)
        self.xsecF.SetNpx(10000)
        
        self.dxsecF = TF1("PositrionAnnihilationDcosT",self.getdXsection,-1,1)
        self.dxsecF.SetNpx(10000)
        
        
#        print mbcolors.OKGREEN,"PositronAnnihilation width of A' for mA=",Ma," mChi=",Mchi,"AlphaD=",alphaD,"is: ",self.width,mbcolors.ENDC
    
    def getApWidth(self):
        return self.width
    #See A.C. notes on 25/3/2019
    
    #at fixed s, the angular distribution (z is the cosine of the angle wrt the beam axis
    def getdXsection(self,z):
        val=0;
        if (self.sdXsec/4-self.Mchi*self.Mchi < 0):
            val=0
        else:
            val=(self.sdXsec/2.+(self.sdXsec/4.-self.Mchi*self.Mchi)*(z[0]*z[0]-1))
        return val
        
    #total xsection vs s
    def getXsectionS(self,s):
        val=0;
        if (s/4-self.Mchi*self.Mchi < 0):
            val=0.
        else:
            val = 1;
            val=val*4*math.pi*PositronAnnihilationKin.alpha*self.alphaD*self.eps*self.eps;
            val=val/math.sqrt(s);
            val=val*math.sqrt(s/4-self.Mchi*self.Mchi);
            val=val/((s-self.Ma*self.Ma)*(s-self.Ma*self.Ma)+s*s*self.width*self.width/(self.Ma*self.Ma))
            val=val*(2.*s/3.+4.*self.Mchi*self.Mchi/3.) #This is in GeV^-2
            val=val*392E6; #this is in pbarn
        return val

    # Here E is the positron energy in lab frame    
    def getXsectionE(self,E):
      #  sValue=2*E[0]*PositronAnnihilationKin.Me+PositronAnnihilationKin.Me*PositronAnnihilationKin.Me
        sValue=2*E[0]*PositronAnnihilationKin.Me;
        return self.getXsectionS(sValue)
    
    def integrateXsectionE(self,Emin,Emax):
        mEmin=Emin
        mEmax=Emax
        if (Emin<self.EA):
            print mbcolors.WARNING,"PositronAnnihilation integrateXsectionE lower bound",Emin,"is less than the value in the initialization",EA,mbcolors.ENDC
            mEmin=self.EA
        if (Emax>self.EB):
            print mbcolors.WARNING,"PositronAnnihilation integrateXsectionE upper bound",Emax,"is greater than the value in the initialization",EB,mbcolors.ENDC
            mEmax=self.EB
        
        return self.xsecF.Integral(mEmin,mEmax);
    
    def integrateXsectionDefRange(self):
        mEmin=self.EA
        mEmax=self.EB
        return self.xsecF.Integral(mEmin,mEmax);
    
    def getRandomCosAngle(self,E):
        self.sdXsec=2*E*PositronAnnihilationKin.Me
        return self.dxsecF.GetRandom();
    
    def getRandomEnergy(self):
        return self.xsecF.GetRandom();
    
class PositronAnnihilationSpectra:
  
    def __init__(self,Ma,Mchi,eps,alphaD,hE):
            self.hE=hE;
            self.Ma=Ma
            self.Mchi=Mchi
            self.eps=eps
            self.alphaD=alphaD
            self.nbins=hE.GetNbinsX()
            self.posKinCalculators = [];   
            self.xsection = [];            #this is Int(bin) sigma(E) dE  ,
            self.xsectionAvg = [];         #this is Int(bin) sigma(E) T(E) dE =  [ Int(bin)sigma(E) ] * T(E) 
            self.mrand=TRandom3(0) 
             
            print mbcolors.OKBLUE,"PositronAnnihilationSpectra there are: ",self.nbins," bins ",mbcolors.ENDC
            for ii in range(self.nbins):
                Ei = hE.GetBinCenter(ii+1);
                deltaE = hE.GetXaxis().GetBinWidth(ii+1) #this is the bin width
                Eimin = Ei-deltaE/2;
                Eimax = Ei+deltaE/2;
                self.posKinCalculators.append(PositronAnnihilationKin(Ma,Mchi,eps,alphaD,Eimin,Eimax))
            self.computeXsections()
    
    
    def computeXsections(self):
        self.xsection=[];
        self.xsectionAvg = []; 
        for ii in range(self.nbins):
            Ei = self.hE.GetBinCenter(ii+1);
            deltaE = self.hE.GetXaxis().GetBinWidth(ii+1) #this is the bin width
            TE = self.hE.GetBinContent(self.hE.FindBin(Ei))
            self.xsection.append(self.posKinCalculators[ii].integrateXsectionDefRange())
            self.xsectionAvg.append(self.xsection[ii]*TE)
            print mbcolors.OKGREEN,"PositronAnnihilationSpectra bin: ",ii,"xsec(pb): ",self.xsection[ii],self.xsectionAvg[ii],TE,deltaE,mbcolors.ENDC   
    
    #This assumes that the current LHE file can be written!  
    #Weight is the weight of each event            
    def generateAndWriteEvents(self,nEvt,ibin,weight,LHEfname):
        kinCalc=self.posKinCalculators[ibin];
        Ei = self.hE.GetBinCenter(ibin+1);
        dest_lhefile=open(LHEfname,"a");
        for iv in range(nEvt):
            Ebeam=(kinCalc.getRandomEnergy());
            z=kinCalc.getRandomCosAngle(Ebeam);
            phi=self.mrand.Uniform(-math.pi,math.pi);       
            s=Ebeam*2*PositronAnnihilationKin.Me
            E=math.sqrt(s)/2.
            q=math.sqrt(E*E-self.Mchi*self.Mchi)
            Pep=TLorentzVector(0,0,math.sqrt(Ebeam*Ebeam-PositronAnnihilationKin.Me*PositronAnnihilationKin.Me),Ebeam)
            Pem=TLorentzVector(0,0,0,math.sqrt(Ebeam*Ebeam-PositronAnnihilationKin.Me*PositronAnnihilationKin.Me),Ebeam)
            
            Pboost=(Pep+Pem).BoostVector()
            
            #chi and chibar in CM
            Pchi=TLorentzVector(q*math.sqrt(1-z*z)*math.sin(phi),q*math.sqrt(1-z*z)*math.cos(phi),q*z,E)
            Pchibar=TLorentzVector(-q*math.sqrt(1-z*z)*math.sin(phi),-q*math.sqrt(1-z*z)*math.cos(phi),-q*z,E)
            
            Pchi.Boost(Pboost)
            Pchibar.Boost(Pboost)        
                    
        dest_lhefile.close()
        
    #This assumes that the current LHE file can be written!  
    #Weight is the weight of each event            
    def generateAndWriteEventsWRotation(self,nEvt,ibin,weight,hAngleTMP,LHEfname):
        kinCalc=self.posKinCalculators[ibin];
        Ei = self.hE.GetBinCenter(ibin+1);
        dest_lhefile=open(LHEfname,"a");
        for iv in range(nEvt):
            Ebeam=(kinCalc.getRandomEnergy());
            z=kinCalc.getRandomCosAngle(Ebeam);
           
            s=Ebeam*2*PositronAnnihilationKin.Me
            E=math.sqrt(s)/2.
            q=math.sqrt(E*E-self.Mchi*self.Mchi)
            
                    
            Pep=TLorentzVector(0,0,math.sqrt(Ebeam*Ebeam-PositronAnnihilationKin.Me*PositronAnnihilationKin.Me),Ebeam)
            Pem=TLorentzVector(0,0,0,PositronAnnihilationKin.Me)
            
            Pboost=(Pep+Pem).BoostVector()
            
            #chi and chibar in CM
            phi=0.
            Pchi=TLorentzVector(q*math.sqrt(1-z*z)*math.sin(phi),q*math.sqrt(1-z*z)*math.cos(phi),q*z,E)
            Pchibar=TLorentzVector(-q*math.sqrt(1-z*z)*math.sin(phi),-q*math.sqrt(1-z*z)*math.cos(phi),-q*z,E)
            
            Pchi.Boost(Pboost)
            Pchibar.Boost(Pboost)
            
            ctheta=hAngleTMP.GetRandom()
            theta=TMath.ACos(ctheta)           #This returns 0 if ctheta > 1
            phi=self.mrand.Uniform(-math.pi,math.pi);      
            Pep.RotateY(theta)
            Pep.RotateZ(phi)
            Pem.RotateY(theta)
            Pem.RotateZ(phi)
            Pchi.RotateY(theta)
            Pchi.RotateZ(phi)
            Pchibar.RotateY(theta)
            Pchibar.RotateZ(phi)
                  
                  
            self.writeLHEevent(dest_lhefile,Pep,Pem,Pchi,Pchibar,nEvt,weight)
        dest_lhefile.close()    
        
    def writeLHEevent(self,dest_lhefile,Pep,Pem,Pchi,Pchibar,nEvt,weight):
        newline="<event>\n"
        dest_lhefile.write(newline)
        #First line: number of particles 0 weightPerParticle totWeight alphaEM alphaS 
        #Tot weight is WRONG, but we don't use it here.
        newline="4 0 "+str(1.*weight)+" "+str(weight*nEvt)+" "+str(1./137)+" "+str(1.18E-1)+"\n"
        dest_lhefile.write(newline)
            
        #Second line: beam e+ 
        #11 -1 0 0 0 0 Px Py Pz E M 0. -1.
        newline="-11 -1 0 0 0 0 "+str(Pep.Px())+" "+str(Pep.Py())+" "+str(Pep.Pz())+" "+str(Pep.E())+" "+str(Pep.M())+" 0. -1."+"\n"
        dest_lhefile.write(newline)
        
        #Third line: e- 
        newline="11 -1 0 0 0 0 "+str(Pem.Px())+" "+str(Pem.Py())+" "+str(Pem.Pz())+" "+str(Pem.E())+" "+str(Pem.M())+" 0. -1."+"\n"
        dest_lhefile.write(newline)
            
        #Fourth line: chi
        #611 1 1 1 0 0 Px Py Pz E M 0. -1.
        newline="611 1 1 1 0 0 "+str(Pchi.Px())+" "+str(Pchi.Py())+" "+str(Pchi.Pz())+" "+str(Pchi.E())+" "+str(Pchi.M())+" 0. -1."+"\n"
        dest_lhefile.write(newline)    
        #Fifth line: chiBar
        #-611 1 1 1 0 0 Px Py Pz E M 0. -1.
        newline="-611 1 1 1 0 0 "+str(Pchibar.Px())+" "+str(Pchibar.Py())+" "+str(Pchibar.Pz())+" "+str(Pchibar.E())+" "+str(Pchibar.M())+" 0. -1."+"\n"   
        dest_lhefile.write(newline)

        newline="</event>\n"
        dest_lhefile.write(newline)
        
        
            
     
     
     
     
                