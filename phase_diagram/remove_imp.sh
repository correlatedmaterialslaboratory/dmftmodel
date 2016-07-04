#!/bin/bash
x=0
while [ $x -le 5 ]
do
	y=0
	while [ $y -le 9 ]
	do
	if [ -d "$x.$y" ]; then
		rm -r "$x.$y"/nohup*
		rm -r "$x.$y"/status*
	fi
	y=$(($y+1))
	done
x=$(($x+1))
done
