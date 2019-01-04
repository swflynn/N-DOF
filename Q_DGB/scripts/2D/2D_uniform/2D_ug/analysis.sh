# run the fortran code
./a.out < input.dat
# generate data analysis
python analysis.py eigenvalues.dat theory.dat
# move everything to new directory for plotting
mkdir data
mv all.dat centers.dat eigenvalues.dat overlap_eigenvalues.dat simulation.dat theory.dat data
cd data
# truncate data files for convenience
head -21 all.dat > plot.dat
head -100 theory.dat > temp.dat
mv temp.dat theory.dat
