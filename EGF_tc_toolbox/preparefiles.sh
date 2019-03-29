#!/bin/bash

export SAC_DISPLAY_COPYRIGHT=0

echo 'Converts between sac and mseed and merges the files into daily files for the given stations and dates, the original files will be permantly changed, the files and script should therefore be copied to a new folder before executing. Input filename must contain stationname $stat and date $dd in spesified dateformat. Output sac file is $network.$stat.00.HH$ch.D.$dd.000000.SAC'

echo 'Spesify the input file format (sac or mseed)'
read fformat

echo 'Spesify network'
read network

echo 'Spesify the stations'
read -a stats

#echo spesify the channels
#read -a channels
ch=Z

echo 'Give the first and last day and dateformat (default is %Y.%j'
read firstday lastday dformat
dformat=${format:-%Y.%j}

echo 'Decimate (give two numbers if the new samplingrate >7, press enter to not decimate)'
read dec1 dec2
dec1=${dec1:-0}
dec2=${dec2:-1}

fd=$(date --date=$firstday +%Y%m%d)
ld=$(date --date=$lastday +%Y%m%d)

for stat in ${stats[@]}
do

#for ch in ${channels[@]}
#do

dy=$fd
while [[ "$dy" -le $ld ]]
do
	
dd=$(date --date=$dy +%Y.%j)
echo $dd
dd2=$(date --date=$dy +$dformat)

if [ "$fformat" == "mseed" ] 
then
mseed2sac $network*$stat*$dd2*
fi

ns=0
for file in *$stat*$dd*SAC
do
# Count the number of files per day
ns=$(($ns+1))

if [ $dec1 -gt 1 ]
then

echo $dec1 $dec2

sac <<EOF
r $file
decimate $dec1 filter on
decimate $dec2 filter on
listhdr delta
write over
chnhdr delta $ndel
listhdr npts
listhdr delta
write $stat-$ch.$year-$mnd-$i.sac-$ns
quit
EOF

else

mv $file $stat-$ch.$dd.sac-$ns

fi

done

echo $ns

if [ $ns -gt 1 ] 
then
 
sac <<EOF
m mergeloop.m number $ns file $stat-$ch.$dd.sac 
quit
EOF

mv $stat-$ch.$dd.sac $stat.00.HH$ch.D.$dd.000000.SAC

rm $stat-$ch.$dd.sac-*

else

mv $stat-$ch.$dd.sac-$ns $network.$stat.00.HH$ch.D.$dd.000000.SAC

fi

dy=$(date -d"$dy + 1 day" +%Y%m%d)

done
done
