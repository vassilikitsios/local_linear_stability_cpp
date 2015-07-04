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
	 	 
# ------------------ Spatial Growth Rate verses Real Frequency Map ----------------------

set xlabel "or"
set xrange[*:*]

set ylabel "-ai" 
set yrange[*:*]

plot '../results/eigen_mapping_sorted.mode_0008.real_omega.dat' using 4:3 with lines lt 1 lw 2 title "",\
  	'../results/eigen_mapping_sorted.mode_0008.real_omega_raw.dat' using 4:3 with points pt 4 title ""


# -------------------------------- Common coloured settings  -------------------------------- 

reset
set size ratio 1.0
set terminal postscript enhanced color font "Times-Roman" 20 dashlength 3
unset key

set style line 100 lt 1 lw 0.5 lc 0
set style line 200 lt 2 lw 2 lc 0

#unset hidden3d 
#unset surf
set pm3d map
set pm3d at b hidden3d 100 
set pm3d corners2color c1

#set palette gray positive
set palette rgbformulae 33,13,10

#set surface
#set contour base
#set cntrparam levels discrete 0.0
#set nokey
#set view 0,0,1,1

unset zlabel
unset ztics

# -------------------------------- Wave Number Map -------------------------------- 

unset label
unset arrow

set xlabel "ar"
#set xrange[1.1:6.1]
#set xtics 2,1,6

set ylabel "ai"
#set yrange[-5.6:0.5]
#set ytics -5,1,0

set cblabel ""
#set label "oi" at 6.23,1.05
#set cbrange [-0.5:0.5]
#set cbtics -0.5,0.25,0.5

splot '../results/eigen_mapping_sorted.mode_0008.dat' using 1:2:6 every :::0::20 with lines ls 100 title "",\
      '../results/eigen_mapping_sorted.mode_0008.real_omega_raw.dat' using 1:2:2 ls 200 title ""

# -------------------------------- Complex Frequency Map --------------------------------

unset label

set xlabel "or"
#set xrange[0:4]
#set xtics 0,1,4

set ylabel "oi"
#set yrange[-1.4:0.8]
#set ytics -1.2,0.4,0.8

set cblabel ""
#set label "ai" at 4.1,1.0
#set cbrange [-5:0]
#set cbtics -5,1,0

#set label "M1" at 1.1,0.6
#set label "M2" at 0.8,-1
#set label "P" at 3.5,-0.8
#set arrow from 3.4,-0.7 to 2.8,-0.45 filled

splot '../results/eigen_mapping_sorted.mode_0008.dat' using 5:6:2 every :::0::20 with lines ls 100 title "k imag"

