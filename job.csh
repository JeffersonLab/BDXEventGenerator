#!/bin/tcsh -f

#$1: tag name.
#File name with the parameters is: run/$1.txt
#TAG mA mCHI eps alphaD Nevents Ebeam Ldump Lx Ly Lz procid
#0   1  2   l 3   4      5       6     7     8  9  10 11


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
    rm -rf $workdir
endif
mkdir $workdir
echo "ls again:"
ls -ltrah
echo "done"

#1: copy everything here in the work dir
cp -Rp $BDX_EVENT_GENERATOR/* $workdir
#2: clean some space, we do not need events 
cd $workdir
mv $workdir/AprimeAlAlpha1/Events/banner_header.txt ./tmp_header
rm -rf $workdir/AprimeAlAlpha1/Events/*
mv ./tmp_header $workdir/AprimeAlAlpha1/Events/banner_header.txt 
#3: prepare the cards
cd $workdir
set fname=run/$1.txt
python -c "import utils;utils.createCards('$fname')"
#echo $?
#4: copy the cards
cp $workdir/run/run_card.dat $workdir/AprimeAlAlpha1/Cards
cp $workdir/run/param_card.dat $workdir/AprimeAlAlpha1/Cards
#5: go AprimeAlAlpha1
cd $workdir/AprimeAlAlpha1
#6 run
./bin/generate_events 0 $1
#7 copy the events we got back to the main place
cd $workdir/AprimeAlAlpha1/Events
cp -rp * $BDX_EVENT_GENERATOR/AprimeAlAlpha1/Events
#8 generate the html
cd $BDX_EVENT_GENERATOR/AprimeAlAlpha1
./bin/gen_crossxhtml-pl $1

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
