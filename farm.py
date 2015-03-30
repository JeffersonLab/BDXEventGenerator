#!/usr/local/bin/python

import string,math,os,sys,subprocess,shutil,utils

EventGeneratorLocation=os.environ['BDX_EVENT_GENERATOR']
os.chdir(EventGeneratorLocation)
runList=open("run/runs.dat", "r")
lines = runList.readlines()
#each line must be written as:
#TAG Ma Mchi epsilon alphaD Nevents Ebeam Ldump Lx Ly Lz
#each line starting with "#" is a comment
command = "rm run/log*"
os.system(command)
print "Doing runs"
print "tag \t mA \t mChi \t eps \t alphaD \t Nevents \t Ebeam \t Ldump \t Lx \t Ly \t Lz \t procid"
for line in lines:
	if ((line[0]!="#")and(len(line)>0)):
		x = line.split()
		if (len(x)!=12):
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
		fname = "run/"+tag+".txt"
		runfile = open(fname,"w")
		runfile.write(x[0]+" "+x[1]+" "+x[2]+" "+x[3]+" "+x[4]+" "+x[5]+" "+x[6]+" "+x[7]+" "+x[8]+" "+x[9]+" "+x[10]+" "+x[11]+"\n")
		runfile.close()
		command = "bsub -q long -P sl5_64 -o run/log"+tag+".o -e run/log"+tag+".e -R \"rusage[mem=500,swp=1000]\" ./job.csh "+tag
		print command
		os.system(command)
