for x in 1 2 3 4 5
do
	f="round$x"
	if [ -d $f ]; then
		continue
	else
		mkdir $f
		cp data $f
		cp Delta_site.dat $f
		cp sub* $f
		mv convergence $f
		break
	fi
done
