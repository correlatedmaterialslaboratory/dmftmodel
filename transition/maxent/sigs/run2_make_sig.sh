N=250
L=10.
mesh=5000

mkdir sig_i
cd sig_i
cp ../Sig_i.out ./
cp ../spline.py ./
python spline.py Sig_i.out Sig_reg.out
cp ../maxent_params.dat ./
cp ../remove_maxent.sh ./
maxent_run.py Sig_reg.out
sed '1,2 d' Sig.out > real_sig_me.out
cd ..

mkdir sig_m
cd sig_m
cp ../Sig_m.out ./
cp ../spline.py ./
python spline.py Sig_m.out Sig_reg.out
cp ../maxent_params.dat ./
cp ../remove_maxent.sh ./
maxent_run.py Sig_reg.out
sed '1,2 d' Sig.out > real_sig_me.out
cd ..
