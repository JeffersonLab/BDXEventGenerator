#!/bin/tcsh -f


#$1: tag name.
#File name with the parameters is: run/$1.txt

if ($#argv != 1) then
        echo "Usage: $0 tag name"
        goto done
endif

set rootfile="2p2GeV.Al.root"
#set rootfile="6GeV.Al.root"
#set rootfile="11GeV.Al.root"
#set rootfile="4p315GeV.Al.root"

echo "go to work dir"
cd /project/Gruppo3/fiber5/celentano/EventGenerator/runDir
echo "ls:"
ls -ltrah
echo "done"
du -sh
echo "done du -sh"
set workdir="`pwd`/EventGenerator_$1"
echo $workdir
if (-e $workdir) then
    echo "remove  already-existing workdir"
    rm -rf $workdir
endif
mkdir $workdir
echo "ls again:"
ls -ltrah
echo "where is ld?:"
which ld
echo "done"

#0: set new variables
setenv BDX_EVENT_GENERATOR_ORIG $BDX_EVENT_GENERATOR
setenv BDX_EVENT_GENERATOR $workdir

echo "new variables"
echo "orig: $BDX_EVENT_GENERATOR_ORIG"
echo "new : $BDX_EVENT_GENERATOR"

#1: copy everything here in the work dir
rsync -rl --exclude .git --exclude run --exclude DetectorInteraction/Events/ --exclude AprimeAlAlpha1/Events/ --exclude html --exclude runDir --exclude log $BDX_EVENT_GENERATOR_ORIG/ $workdir
#2: missing dirs
cd $workdir/AprimeAlAlpha1
mkdir Events
cd Events
cp $BDX_EVENT_GENERATOR_ORIG/AprimeAlAlpha1/Events/banner_header.txt ./
cd $workdir/DetectorInteraction
mkdir Events

#3: prepare the cards
cd $workdir
mkdir run
set fname=$BDX_EVENT_GENERATOR_ORIG/run/$1.txt
python -c "import CardsUtils;CardsUtils.createCards('$fname')"
#echo $?
#4: copy the cards
cp $workdir/run/run_card.dat $workdir/Cards
cp $workdir/run/param_card.dat $workdir/Cards
#5: go workdir
cd $workdir
#6 run
echo "run it now"
echo "this is the script"
cat RunEventGenerator.py
echo "this is the run card"
cat $workdir/Cards/run_card.dat
echo "this is the param card"
cat $workdir/Cards/param_card.dat
echo "RUN RUN RUN"
python -u RunEventGenerator.py --run_name $1 --run_card $workdir/Cards/run_card.dat --param_card $workdir/Cards/param_card.dat --root_file $rootfile 
#7 copy the events we got back to the main place
cd $workdir/DetectorInteraction/Events
cp -rp * $BDX_EVENT_GENERATOR_ORIG/DetectorInteraction/Events

#9 remove these temporary location
echo "CLEAR"
cd $workdir
cd ..
ls
rm -rf $workdir
ls
goto done

done:
    exit(0)
