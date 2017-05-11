reset
set terminal eps
set output "plot_n_orbs.eps"
plot "<awk '$2==0.0002 && $4==2e-6 && $3==3582' out_1126685.dat" using (1/$1):6
