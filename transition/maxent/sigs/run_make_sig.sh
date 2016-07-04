N=250
L=10.
mesh=5000
for x in {9..19}
do
	mkdir sig.$x
	cd sig.$x
	cp ../Sig.out.$x ./
	cp ../spline.py ./
	python spline.py Sig.out.$x Sig_reg.out
	#pade_sig.py Sig_reg.out -L $L -N $N -mesh $mesh -s real_sig.out
	cp ../maxent_params.dat ./
	cp ../remove_maxent.sh ./
	maxent_run.py Sig_reg.out
	sed '1,2 d' Sig.out > real_sig_me.out
	cd ..
done
