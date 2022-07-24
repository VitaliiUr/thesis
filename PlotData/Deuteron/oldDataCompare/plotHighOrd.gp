set term pdf enhanced





  set style line 1 lt rgb "orange" lw 3 dt 1
  set style line 2 lt rgb "yellow" lw 3 dt 4
  set style line 3 lc rgb "dark-green" lw 3 dt 2
  set style line 4 lc rgb "dark-violet" lw 3 dt 3
  set style line 5 lc rgb "blue" lw 3 dt 5
  # set style line 6 lc rgb "violet" lw 1.5 dt 5
  # set style line 7 lc rgb "cyan" lw 1.5  pt 3 pi -8 ps 0.3
  set pointintervalbox 0.02  ## interval to a point




set output "HighOrd100MeV.pdf"

set tics scale 2
set mxtics 4
set mytics 2


set ylabel "d^2{/Symbol s}/d{/Symbol W} [{/Symbol m}b sr^{-1}]" offset 1.3,0 
set xlabel "{/Symbol q}_p [deg]"	
set key Left at 175, graph 0.95  
# unset key

plot "../data/CROSS_100mev.dat" u 1:5 w l ls 3 title "SMS N3LO",\
 "../data/CROSS_100mev.dat" u 1:7 w l ls 4 title "SMS N4LO",\
 "../data/CROSS_100mev.dat" u 1:10 w l  ls 5 title "SMS N4LO+",\
  "../dataSCS/CROSS_old_100mev.dat" u 1:5 w l ls 1 title "SCS N3LO",\
  "../dataSCS/CROSS_old_100mev.dat" u 1:6 w l ls 2 lc "red" title "SCS N4LO",\
  "../ExpData/Exp100.dat" u 1:2:3 notitle pt 7 ps 0.5 lc rgb "black" with yerrorbars,\
  "../ExpData/Exp100_2.dat" u 1:2:3 notitle pt 6 ps 0.5 lc rgb "black" with yerrorbars,\
  "../ExpData/Exp100n.dat" u 1:2:3 notitle pt 5 ps 0.5 lc rgb "black" with yerrorbars,\
  "../ExpData/Exp100_3.dat" u 1:2:3 notitle pt 9 ps 0.5 lc rgb "black" with yerrorbars




set output "HighOrd30MeV.pdf"


set multiplot

set ylabel "d^2{/Symbol s}/d{/Symbol W} [{/Symbol m}b sr^{-1}]"
set xlabel "{/Symbol q}_p [deg]"	


plot "../data/CROSS_30mev.dat" u 1:5 w l ls 3 title "SMS N3LO",\
 "../data/CROSS_30mev.dat" u 1:7 w l ls 4 title "SMS N4LO",\
 "../data/CROSS_30mev.dat" u 1:10 w l  ls 5 title "SMS N4LO+",\
  "../dataSCS/CROSS_old_30mev.dat" u 1:5 w l ls 1 title "SCS N3LO",\
  "../dataSCS/CROSS_old_30mev.dat" u 1:6 w l ls 2 lc 'red' title "SCS N4LO",\
  "../ExpData/Exp30.dat" u 1:2:3 notitle pt 7 ps 0.5 lc rgb "black" with yerrorbars,\
   "../ExpData/Exp30_2.dat" u 1:2:3 notitle pt 6 ps 0.5 lc rgb "black" with yerrorbars


set nokey
unset ylabel
unset xlabel
set origin 0.3, 0.2
set size 0.35, 0.45
set xrange [70:90]
set yrange [34:38]
set xtics 10
set ytics 1
set mytics 2
plot "../data/CROSS_30mev.dat" u 1:5 w l ls 3 title "SMS N3LO",\
 "../data/CROSS_30mev.dat" u 1:7 w l ls 4 title "SMS N4LO",\
 "../data/CROSS_30mev.dat" u 1:10 w l  ls 5 title "SMS N4LO+",\
  "../dataSCS/CROSS_old_30mev.dat" u 1:5 w l ls 1 title "SCS N3LO",\
  "../dataSCS/CROSS_old_30mev.dat" u 1:6 w l ls 2 lc 'red' title "SCS N4LO"

unset multiplot



# set nokey
# unset ylabel
# unset xlabel
# set origin 0.08, 0.1
# set size 0.35, 0.45
# set xrange [55:65]
# set yrange [5.5:6.5]
# set xtics 5
# set ytics 0.2
# set mytics 2
# set mxtics 2
# replot
