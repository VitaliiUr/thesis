set term pdf enhanced


MP_LEFT = .1
MP_RIGHT = .95
MP_BOTTOM = .17
MP_TOP = .95
MP_GAP = 0.0005


  set style line 1 lt rgb "orange" lw 3 dt 1
  set style line 2 lt rgb "red" lw 3 dt 4
  set style line 3 lc rgb "dark-green" lw 3 dt 2
  set style line 4 lc rgb "dark-violet" lw 3 dt 3
  set style line 5 lc rgb "blue" lw 3 dt 5
  # set style line 6 lc rgb "violet" lw 1.5 dt 5
  # set style line 7 lc rgb "cyan" lw 1.5  pt 3 pi -8 ps 0.3
  set pointintervalbox 0.02  ## interval to a point




set output "HighOrd100MeV_trunc.pdf"

# set multiplot layout 1,2 \
#               margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP
set multiplot layout 1,2 \
              margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP

set tics scale 2
set mxtics 4
set mytics 2
set xtics 0,20,179

set ylabel "d^2{/Symbol s}/d{/Symbol W} [{/Symbol m}b sr^{-1}]" offset 1.3,0 
set xlabel "{/Symbol q}_p [deg]"	
# set key Left at 175, graph 0.95  
set key Left at 110, graph 0.35
# unset key
set label 1 "(a)" at 135, graph 0.8 font ",14"

plot "../data/CROSS_100mev.dat" u 1:5 w l ls 3 title "SMS N3LO",\
 "../data/CROSS_100mev.dat" u 1:7 w l ls 4 title "SMS N4LO",\
 "../data/CROSS_100mev.dat" u 1:10 w l  ls 5 title "SMS N4LO+",\
  "../dataSCS/CROSS_old_100mev.dat" u 1:5 w l ls 1 title "SCS N3LO",\
  "../dataSCS/CROSS_old_100mev.dat" u 1:6 w l ls 2 title "SCS N4LO",\
  "../ExpData/Exp100.dat" u 1:2:3 notitle pt 7 ps 0.5 lc rgb "black" with yerrorbars,\
  "../ExpData/Exp100_2.dat" u 1:2:3 notitle pt 6 ps 0.5 lc rgb "black" with yerrorbars,\
  "../ExpData/Exp100n.dat" u 1:2:3 notitle pt 5 ps 0.5 lc rgb "black" with yerrorbars,\
  "../ExpData/Exp100_3.dat" u 1:2:3 notitle pt 9 ps 0.5 lc rgb "black" with yerrorbars

unset ylabel
set yrange [1:8]
set ytics format ""
set xtics 0,20,180
set key Left at 110, graph 0.3

set label 1 "(b)"

set style fill transparent solid 0.25 # partial transparency
# set style fill noborder # no separate top/bottom lines
plot "../data/Trunc_VU100MeV_cross_sms.dat" u 1:6  notitle "SMS N^4LO+" lc "blue" lw 0.0001 w filledcurves,\
      "../data/CROSS_100mev.dat" u 1:10 w l ls 5 lw 1  title "SMS N4LO+",\
      "../dataSCS/Trunc_VU100MeV_Cross_scs.dat" u 1:5 notitle "SCS N^4LO" lw 0.0001  lc "red"    w filledcurves,\
      "../dataSCS/CROSS_old_100mev.dat" u 1:6 w l  ls 2 lw 1 lc "red"  title "SCS N4LO"

unset multiplot


set output "HighOrd30MeV_trunc.pdf"


set multiplot layout 1,2
              # margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP
# MP_LEFT = .1
# MP_RIGHT = .95
# MP_BOTTOM = .17
# MP_TOP = .95
# MP_GAP = 0.0005
set origin 0.0, 0.0
set size 0.53, 1.0
set lmargin 6
set rmargin 0
set tmargin 1
set bmargin 3
set yrange [0:40]
set xtics 0,20,179
set nokey

set ylabel "d^2{/Symbol s}/d{/Symbol W} [{/Symbol m}b sr^{-1}]"
set xlabel "{/Symbol q}_p [deg]" offset 0,0.5	
set format y "%1.f"
set label 1 "(a)" at 135, graph 0.8 font ",14" 

plot "../data/CROSS_30mev.dat" u 1:5 w l ls 3 title "SMS N3LO",\
 "../data/CROSS_30mev.dat" u 1:7 w l ls 4 title "SMS N4LO",\
 "../data/CROSS_30mev.dat" u 1:10 w l  ls 5 title "SMS N4LO+",\
  "../dataSCS/CROSS_old_30mev.dat" u 1:5 w l ls 1 title "SCS N3LO",\
  "../dataSCS/CROSS_old_30mev.dat" u 1:6 w l ls 2 title "SCS N4LO",\
  "../ExpData/Exp30.dat" u 1:2:3 notitle pt 7 ps 0.5 lc rgb "black" with yerrorbars,\
   "../ExpData/Exp30_2.dat" u 1:2:3 notitle pt 6 ps 0.5 lc rgb "black" with yerrorbars

set lmargin 0
set rmargin 2
set tmargin 1
set bmargin 3
set autoscale
set origin 0.53, 0.0
set size 0.45, 1.0
unset ylabel
set ytics format ""
set xtics 0,20,180
# set key Left at 110, graph 0.3
set label 1 "(b)" 
set nokey

set style fill transparent solid 1 # partial transparency
# set style fill noborder # no separate top/bottom lines
plot "../data/Trunc_VU30MeV_cross_sms.dat" u 1:6  notitle "SMS N^4LO+" lc "blue" lw 0.0001 w filledcurves,\
      "../data/CROSS_30mev.dat" u 1:10 w l ls 5 lw 1  title "SMS N4LO+",\
      "../dataSCS/Trunc_VU30MeV_Cross_scs.dat" u 1:5 notitle "SCS N^4LO" lw 0.0001  lc "red"  w filledcurves,\
      "../dataSCS/CROSS_old_30mev.dat" u 1:6 w l  ls 2 lw 1 lc "red"  title "SCS N4LO"

unset label 1
unset ylabel
unset xlabel
set format y "%1.f"
set origin 0.213, 0.12  
set size 0.195, 0.46
set xrange [70:90]
set yrange [34:38]
set xtics 10
set ytics 2 offset 0.5,0
set mytics 2
plot "../data/CROSS_30mev.dat" u 1:5 w l ls 3 title "SMS N3LO",\
 "../data/CROSS_30mev.dat" u 1:7 w l ls 4 title "SMS N4LO",\
 "../data/CROSS_30mev.dat" u 1:10 w l  ls 5 title "SMS N4LO+",\
  "../dataSCS/CROSS_old_30mev.dat" u 1:5 w l ls 1 title "SCS N3LO",\
  "../dataSCS/CROSS_old_30mev.dat" u 1:6 w l ls 2 title "SCS N4LO"

set origin 0.655, 0.12  
set size 0.195, 0.46
set xrange [70:90]
set yrange [34:38]
set xtics 10
set ytics 2 offset 0.5,0
set mytics 2
plot "../data/Trunc_VU30MeV_cross_sms.dat" u 1:6  notitle "SMS N^4LO+" lc "blue" lw 0.0001 w filledcurves,\
      "../data/CROSS_30mev.dat" u 1:10 w l ls 5 lw 1  title "SMS N4LO+",\
      "../dataSCS/Trunc_VU30MeV_Cross_scs.dat" u 1:5 notitle "SCS N^4LO" lw 0.0001  lc "red"   w filledcurves,\
      "../dataSCS/CROSS_old_30mev.dat" u 1:6 w l  ls 2 lw 1 lc "red" title "SCS N4LO"



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
