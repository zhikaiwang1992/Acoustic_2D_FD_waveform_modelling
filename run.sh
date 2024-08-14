#! /bin/sh
 NZ=301
 NX=301
 ABC=40
 Pout=../outdata/
 Vin=../inputdata/VHome2000.bin
./main nx=$NX nz=$NZ nt=401 sx=151 sz=151 dx=10.0 dz=10.0 dt=0.001 f=8 flag=16 abc=$ABC out=$Pout vel=$Vin
#ximage n1=$NZ  n2=$NX perc=99 wbox=$NX hbox=$NZ<  $Pout &
