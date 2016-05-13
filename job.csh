#!/bin/tcsh -f

#$1: tag name.
#File name with the parameters is: run/$1.txt
#TAG mA mCHI eps alphaD Nevents Ebeam Ldump Lx Ly Lz procid
#0   1  2     3   4      5       6     7     8  9  10 11


if ($#argv != 1) then
        echo "Usage: $0 tag name"
        goto done
endif

echo "go to scrach"
cd /scratch
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
echo "done"

#0: set new variables
setenv BDX_EVENT_GENERATOR_ORIG $BDX_EVENT_GENERATOR
setenv BDX_EVENT_GENERATOR $workdir

echo "new variables"
echo "orig: $BDX_EVENT_GENERATOR_ORIG"
echo "new : $BDX_EVENT_GENERATOR"

#1: copy everything here in the work dir
cp -Rp $BDX_EVENT_GENERATOR_ORIG/* $workdir
#2: clean some space, we do not need events 
cd $workdir
mv $workdir/AprimeAlAlpha1/Events/banner_header.txt ./tmp_header
rm -rf $workdir/AprimeAlAlpha1/Events/*
mv ./tmp_header $workdir/AprimeAlAlpha1/Events/banner_header.txt 
rm -rf $workdir/DetectorInteraction/Events/*
echo "ls again"
ls $workdir
ls $workdir/AprimeAlAlpha1

#3: prepare the cards
cd $workdir
set fname=run/$1.txt
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
python -u RunEventGenerator.py --run_name $1 --run_card $workdir/Cards/run_card.dat --param_card $workdir/Cards/param_card.dat
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
