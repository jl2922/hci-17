set output "plot_eps.eps"
set terminal postscript
set multiplot layout 1,2 \
  title "Correlation Energy vs Epsilons \
  (Slice at 1030 variational dets and 3582 perturbation dets)"

# Left plot.
set xlabel "{/Symbol e}_{var}"
set ylabel "correlation energy (eV)"
set xrange [0:0.001]
set yrange [-0.596:-0.588]
set format x "%.4f"
set format y "%.3f"
set grid

eps_pts = "2.000e-06 4.000e-06 6.000e-06 8.000e-06 1.000e-05"
plot \
  for [i = 1:words(eps_pts)] \
    "<awk '$1==1030 && $3==3582 && $4==".word(eps_pts, i)."' out_1126685.dat" \
    using 2:6 \
    title "{/Symbol e}_{pt} = ".word(eps_pts, i) \
    with points pointtype i, \
  for [i = 1:words(eps_pts)] \
    -5.948156137145e-001  \
    -2.645880272925e-001  / 1030 \
    -7.186688559442e-001  * x \
    +4.254156714498e+000 / 3582 \
    -5.713191123824e+001 * word(eps_pts, i) \
    +4.566804210414e+001 / 1030**2 \
    -5.393204353920e+000 / 1030 * x \
    +6.796702189220e+000 / 1030 / 3582 \
    -2.049793683650e+004 / 1030 * word(eps_pts, i) \
    +4.786844460836e+002 * x**2 \
    +1.491960461547e+001 * x / 3582 \
    -4.255137588561e+004 * x * word(eps_pts, i) \
    +5.094900531321e+002 * (1/3582)**2 \
    +4.502298247703e+004 / 3582 * word(eps_pts, i) \
    +4.419703584713e+006 * (word(eps_pts, i))**2 \
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
    -5.948156137145e-001  \
    -2.645880272925e-001  / 1030 \
    -7.186688559442e-001  * word(eps_vars, i) \
    +4.254156714498e+000 / 3582 \
    -5.713191123824e+001 * x \
    +4.566804210414e+001 / 1030**2 \
    -5.393204353920e+000 / 1030 * word(eps_vars, i) \
    +6.796702189220e+000 / 1030 / 3582 \
    -2.049793683650e+004 / 1030 * x \
    +4.786844460836e+002 * word(eps_vars, i)**2 \
    +1.491960461547e+001 * word(eps_vars, i) / 3582 \
    -4.255137588561e+004 * word(eps_vars, i) * x \
    +5.094900531321e+002 * (1/3582)**2 \
    +4.502298247703e+004 / 3582 * x \
    +4.419703584713e+006 * x**2 \
    notitle \
    linestyle 2
