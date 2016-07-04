mkdir convergence
for x in 20 40 60 80 90 100
do
	python get_do_site.py $(($x-10)) $x > convergence/docc.dat.$x
	python get_nf_site.py $(($x-10)) $x > convergence/nf.dat.$x
done
