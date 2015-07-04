# ---------------------------------------------------------------- 
# GNUPLOT SCRIPT
# usage : gnuplot plot_results.gp
# ---------------------------------------------------------------- 

# -------------------------------- Common settings -------------------------------- 
set terminal postscript enhanced font "Times-Roman" 32 dashlength 3
set title ""
set output "results.eps"

# -------------------------------- Base flow -------------------------------- 

plot '../results/profile.dat' using 3:2 with lines lt 1 title "u"
plot '../results/profile.dat' using 4:2 with lines lt 1 title "u_{y}"
plot '../results/profile.dat' using 5:2 with lines lt 1 title "u_{yy}"

plot '../results/profile.dat' using 13:2 with lines lt 1 title "Ek"
plot '../results/profile.dat' using 14:2 with lines lt 1 title "Ek_{y}"
plot '../results/profile.dat' using 15:2 with lines lt 1 title "Ek_{yy}"

plot '../results/profile.dat' using ($7/$6):2 with lines lt 1 title "nuE1"
plot '../results/profile.dat' using ($9/$6):2 with lines lt 1 title "constE2"
plot '../results/profile.dat' using ($11/$6):2 with lines lt 1 title "constE3"

# -------------------------------- Eigenvalue spectrum -------------------------------- 

set xlabel 'c_r'
set ylabel 'c_i'
set key bottom left
check_error(error) = error > 1.0e-3 ? 1/0 : 1

set xrange[-0.1:1.3]
set yrange[-1.5:0.3]
plot '../results/kx_0001_0001.eigen_values.dat' using ($10/$2):($11/$2*check_error($17)/(1-$15)/$16) with points pt 6 ps 1 title "u S",\
     '../results/kx_0001_0001.eigen_values.dat' using ($10/$2):($11/$2*check_error($17)/(1-$15)/(1-$16)) with points pt 7 ps 0.5 title "u aS",\
     '../results/kx_0001_0001.eigen_values.dat' using ($10/$2):($11/$2*check_error($17)/$15/$16) with points pt 4 ps 1 title "w S",\
     '../results/kx_0001_0001.eigen_values.dat' using ($10/$2):($11/$2*check_error($17)/$15/(1-$16)) with points pt 5 ps 0.5 title "w aS",\
     0 with lines lt 2 title ""
	 	 
# -------------------------------- Eigenvector -------------------------------- 

set xrange[*:*]
set yrange[*:*]

set xlabel "profiles"
set ylabel "wall normal dimension"
plot '../results/kx_0001_0001.eigen_vector.mode_0006.dat' using 2:1 with lines lt 1 title "u_r",\
     '../results/kx_0001_0001.eigen_vector.mode_0006.dat' using 3:1 with lines lt 2 title "u_i"

plot '../results/kx_0001_0001.eigen_vector.mode_0006.dat' using 6:1 with lines lt 1 title "v_r",\
     '../results/kx_0001_0001.eigen_vector.mode_0006.dat' using 7:1 with lines lt 2 title "v_i"

plot '../results/kx_0001_0001.eigen_vector.mode_0006.dat' using 10:1 with lines lt 1 title "w_r",\
     '../results/kx_0001_0001.eigen_vector.mode_0006.dat' using 11:1 with lines lt 2 title "w_i"

plot '../results/kx_0001_0001.eigen_vector.mode_0006.dat' using 14:1 with lines lt 1 title "p_r",\
     '../results/kx_0001_0001.eigen_vector.mode_0006.dat' using 15:1 with lines lt 2 title "p_i"

