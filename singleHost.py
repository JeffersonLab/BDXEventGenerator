#!/usr/local/bin/python

import string,math,os,sys,subprocess,shutil,utils

EventGeneratorLocation=os.environ['BDX_EVENT_GENERATOR']
os.chdir(EventGeneratorLocation)
runList=open("run/runs.dat", "r")
lines = runList.readlines()
#each line must be written as:
#TAG Ma Mchi epsilon alphaD Nevents Ebeam Ldump Lx Ly Lz
#each line starting with "#" is a comment
print "Doing runs"
print "tag \t mA \t mChi \t eps \t alphaD \t Nevents \t Ebeam \t Ldump \t Lx \t Ly \t Lz \t procid"
for line in lines:
	if ((line[0]!="#")and(len(line)>0)):
		x = line.split()
		if (len(x)!=12):
			print("Error with this line: ")
			print(line)
			continue
		tag    = x[0]
		mA     = x[1]
		mChi   = x[2]
		eps    = x[3]
		alphaD = x[4]
		nevents= x[5]
		eBeam  = x[6]
		ldet   = x[7]
		Lx     = x[8]
		Ly     = x[9]
		Lz     = x[10]
		procid = x[11]
		for word in x:
			word=word+"\t"
			sys.stdout.write(word)
		print ""
		utils.createRunCard(nevents,eBeam,ldet,Lx,Ly,Lz,procid)
		utils.createParamCard(mA,mChi,eps,alphaD)
		runfile = open("run/tmp.csh", "w")
		runfile.write("#!/bin/tcsh -f\n")
		runfile.write("cd "+EventGeneratorLocation+" \n")
		runfile.write("cp run/run_card.dat AprimeAlAlpha1/Cards\n")
		runfile.write("cp run/param_card.dat AprimeAlAlpha1/Cards\n")
		runfile.write("cd AprimeAlAlpha1 \n")
		runfile.write("pwd \n")
		runfile.write("./bin/generate_events 0 "+tag+"\n")
		runfile.close()
		os.chmod("run/tmp.csh",0755)
		subprocess.call("run/tmp.csh",shell=True)
		 
