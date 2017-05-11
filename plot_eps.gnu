set output "plot_eps.eps"
set terminal postscript
set multiplot layout 1,2 \
  title "Correlation Energy vs Epsilons \
  (Slice at 1030 variational dets and 3582 perturbation dets)"

# Left plot.
set xlabel "{/Symbol e}_{var}"
set ylabel "correlation energy (eV)"
set xrange [0:0.001]
set yrange [-0.5948:-0.5936]
set format x "%.4f"
set format y "%.4f"
set grid

eps_pts = "2.000e-06 4.000e-06 6.000e-06 8.000e-06 1.000e-05"
plot \
  for [i = 1:words(eps_pts)] \
    "<awk '$1==1030 && $3==3582 && $4==".word(eps_pts, i)."' out_1126685.dat" \
    using 2:6 \
    title "{/Symbol e}_{pt} = ".word(eps_pts, i) \
    with points pointtype i, \
  for [i = 1:words(eps_pts)] \
    -0.5948019807 \
    -0.2696951881 / 1030 \
    -0.7135116468 * x \
    +4.172611948 / 3582 \
    -51.86142507 * word(eps_pts, i) \
    +45.85340983 / 1030**2 \
    -5.103536965 / 1030 * x \
    +10.34682749 / 1030 / 3582 \
    -20151.74143 / 1030 * word(eps_pts, i) \
    +469.7959964 * x**2 \
    +15.47834283 * x / 3582 \
    -41764.78122 * x * word(eps_pts, i) \
    +570.4539941 * (1/3582)**2 \
    +44852.70315 / 3582 * word(eps_pts, i) \
    +3897562.415 * (word(eps_pts, i))**2 \
    notitle \
    linestyle 2

# Right plot.
unset ylabel
set xlabel "{/Symbol e}_{pt}"
set xrange [0:1.0e-5]
set format x "%.0e"
eps_vars = "0.0002 0.0004 0.0006 0.0008 0.0010"
plot \
  for [i = 1:words(eps_vars)] \
    "<awk '$1==1030 && $3==3582 && $2==".word(eps_vars, i)."' out_1126685.dat" \
    using 4:6 \
    title "{/Symbol e}_{var} = ".word(eps_vars, i) \
    with points pointtype i + 5, \
  for [i = 1:words(eps_pts)] \
    -0.5948019807 \
    -0.2696951881 / 1030 \
    -0.7135116468 * word(eps_vars, i) \
    +4.172611948 / 3582 \
    -51.86142507 * x \
    +45.85340983 / 1030**2 \
    -5.103536965 / 1030 * word(eps_vars, i) \
    +10.34682749 / 1030 / 3582 \
    -20151.74143 / 1030 * x \
    +469.7959964 * word(eps_vars, i)**2 \
    +15.47834283 * word(eps_vars, i) / 3582 \
    -41764.78122 * word(eps_vars, i) * x \
    +570.4539941 * (1/3582)**2 \
    +44852.70315 / 3582 * x \
    +3897562.415 * x**2 \
    notitle \
    linestyle 2
