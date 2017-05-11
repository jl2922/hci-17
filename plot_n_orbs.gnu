set output "plot_n_orbs.eps"
set terminal postscript
set multiplot layout 1,2 \
  title "Correlation Energy VS Number of Orbitals \
  (Slice at {/Symbol e}_{var}=0.0002, {/Symbol e}_{pt}=2.0e-6)"

# Left plot.
set xlabel "1 / number of perturbation orbitals"
set ylabel "correlation energy (eV)"
set xrange [0:0.001]
set yrange [-0.5960:-0.5880]
set format x "%.4f"
set format y "%.3f"
set grid
var_orbs = "246 358 514 778 1030"
plot \
  for [i = 1:words(var_orbs)] \
    "<awk '$2==0.0002 && $4==2e-6 && $1==".word(var_orbs, i)."' out_1126685.dat" \
    using (1/$3):6 \
    title word(var_orbs, i)." variational orbitals" \
    with points pointtype i, \
  for [i = 1:words(var_orbs)] \
    -0.5948019807 \
    -0.2696951881 / word(var_orbs, i) \
    -0.7135116468 * 0.0002 \
    +4.172611948 * x \
    -51.86142507 * 2e-6 \
    +45.85340983 / word(var_orbs, i)**2 \
    -5.103536965 / word(var_orbs, i) * 0.0002 \
    +10.34682749 / word(var_orbs, i) * x \
    -20151.74143 / word(var_orbs, i) * 2e-6 \
    +469.7959964 * 0.0002 * 0.0002 \
    +15.47834283 * 0.0002 * x \
    -41764.78122 * 0.0002 * 2e-6 \
    +570.4539941 * x**2 \
    +44852.70315 * x * 2e-6 \
    +3897562.415 * (2e-6)**2 \
    notitle \
    linestyle 2, \
  for [i = 1:words(var_orbs)] \
    "<awk '$2==0.0002 && $4==2e-6 && $1==".word(var_orbs, i)."' out_1126685.dat" \
    using (1/$3):( \
      -0.5948019807 \
      -0.7135116468 * 0.0002 \
      +4.172611948 * (1/$3) \
      -51.86142507 * 2e-6 \
      +469.7959964 * 0.0002 * 0.0002 \
      +15.47834283 * 0.0002 * (1/$3) \
      -41764.78122 * 0.0002 * 2e-6 \
      +570.4539941 * (1/$3)**2 \
      +44852.70315 * (1/$3) * 2e-6 \
      +3897562.415 * (2e-6)**2 \
    ) \
    notitle \
    with points pointtype 12, \
  -0.5948019807 \
    -0.7135116468 * 0.0002 \
    +4.172611948 * x \
    -51.86142507 * 2e-6 \
    +469.7959964 * 0.0002 * 0.0002 \
    +15.47834283 * 0.0002 * x \
    -41764.78122 * 0.0002 * 2e-6 \
    +570.4539941 * x**2 \
    +44852.70315 * x * 2e-6 \
    +3897562.415 * (2e-6)**2 \
    title "Infinite variational orbitals" \
    linestyle 1

# Right plot.
unset ylabel
set xlabel "1 / number of variational orbitals"
set xrange [0:0.005]
pt_orbs = "1030 1478 1850 2378 2838 3582"
plot \
  for [i = 1:words(pt_orbs)] \
    "<awk '$2==0.0002 && $4==2e-6 && $3==".word(pt_orbs, i)."' out_1126685.dat" \
    using (1/$1):6 \
    title word(pt_orbs, i)." perturbation orbitals" \
    with points pointtype i + 5, \
  for [i = 1:words(pt_orbs)] \
    -0.5948019807 \
    -0.2696951881 * x \
    -0.7135116468 * 0.0002 \
    +4.172611948 / word(pt_orbs, i) \
    -51.86142507 * 2e-6 \
    +45.85340983 * x**2 \
    -5.103536965 * x * 0.0002 \
    +10.34682749 * x / word(pt_orbs, i) \
    -20151.74143 * x * 2e-6 \
    +469.7959964 * 0.0002 * 0.0002 \
    +15.47834283 * 0.0002 / word(pt_orbs, i) \
    -41764.78122 * 0.0002 * 2e-6 \
    +570.4539941 / word(pt_orbs, i)**2 \
    +44852.70315 / word(pt_orbs, i) * 2e-6 \
    +3897562.415 * (2e-6)**2 \
    notitle \
    linestyle 2, \
  for [i = 1:words(pt_orbs)] \
    "<awk '$2==0.0002 && $4==2e-6 && $3==".word(pt_orbs, i)."' out_1126685.dat" \
    using (1/$1):( \
      -0.5948019807 \
      -0.2696951881 * (1/$1) \
      -0.7135116468 * 0.0002 \
      -51.86142507 * 2e-6 \
      +45.85340983 * (1/$1)**2 \
      -5.103536965 * (1/$1) * 0.0002 \
      -20151.74143 * (1/$1) * 2e-6 \
      +469.7959964 * 0.0002 * 0.0002 \
      -41764.78122 * 0.0002 * 2e-6 \
      +3897562.415 * (2e-6)**2 \
    ) \
    notitle \
    with points pointtype 12, \
  -0.5948019807 \
    -0.2696951881 * x \
    -0.7135116468 * 0.0002 \
    -51.86142507 * 2e-6 \
    +45.85340983 * x**2 \
    -5.103536965 * x * 0.0002 \
    -20151.74143 * x * 2e-6 \
    +469.7959964 * 0.0002 * 0.0002 \
    -41764.78122 * 0.0002 * 2e-6 \
    +3897562.415 * (2e-6)**2 \
    title "Infinite perturbation orbitals" \
    linestyle 1
#  for [i = 1:words(pt_orbs)] \
#    -5.951369e-01 \
#    -2.406836e-01 * x \
#    -3.349253e-01 * 0.0002 \
#    +5.052496e+00 / word(pt_orbs, i) \
#    -4.256744e+01 * 2e-6 \
#    notitle \
#    linestyle 2, \