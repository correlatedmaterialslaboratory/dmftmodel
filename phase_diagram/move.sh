#take out actual files in the job folder to current folder

x=0
prefix="1805"
fold="naivephya"
while [ $x -le 100 ]
do
	if [ -d "$prefix$x$fold" ]; then
		y=0
		while [ $y -le 6 ]
		do
		mv "$prefix"$x$fold/$y.* ./
		#echo "1448"$x$fold/$y.*
		y=$(($y+1))
		done
	fi
	x=$(($x+1))
	#echo $x
done


