mkdir data1
cd data1
cp ../Sig.out ./
cp ../spline.py ./
python spline.py Sig.out Sig_reg.out
cp ../maxent_params.dat ./
cp ../remove_maxent.sh ./
maxent_run.py Sig_reg.out
sed '1,2 d' Sig.out > real_sig_me.out
cd ..
