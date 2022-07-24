set term pdf enhanced size 5.00in, 2.70in

  set style line 1 lt rgb "orange" lw 2 dt 1 pt 13 pi -10 ps 0.3 ## LO
  set style line 2 lt rgb "blue" lw 2 dt 1  pt 7 pi -10 ps 0.3 ## NLO
  set style line 3 lc rgb "green" lw 2 dt 2 ## N2LO
  set style line 4 lc rgb "red" lw 2 dt 3 ## N3LO
  set style line 5 lc rgb "black" lw 2 dt 1 ## N4LO
  set style line 6 lc rgb "cyan" lw 2 dt 5  ## N4LO+
  set style line 7 lc rgb "violet" lw 2  dt 4 ## AV18
  set pointintervalbox 0.02  ## interval to a point


MP_LEFT = .1
MP_RIGHT = .95
MP_BOTTOM = .17
MP_TOP = .95
MP_GAP = 0.0005

set output 'Py_100MeV_all_trunc.pdf'
set multiplot layout 1,3 \
              margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP

# set size square
#set format y "%.1f"
#set key box opaque
set tics scale 2
set yrange [-0.35:0.12]
set ylabel 'P_y' offset 1.0,0
set xlabel " {/Symbol q}_p [deg]"
set key at 95, graph 0.4 font ",9"
# unset key
set ytics 0.1
set xtics 0,30,179
set mytics 2
set mxtics 3
set label 1 "(a)" at 130, graph 0.8 font ",14"


plot "./data/ProtonPolDeg_100mev.dat" u 1:2  w lp title "LO" ls 1,\
	"" u 1:3  w lp title "NLO" ls 2,\
	"" u 1:4  w l title "N^2LO" ls 3,\
	"" u 1:5  w l title "N^3LO" ls 4,\
	"" u 1:7  w l title "N^4LO" ls 5,\
	"" u 1:10  w l title "N^4LO+" ls 6,\
	 "" u 1:11  w l title "AV18" ls 7

unset ylabel
set format y ""
unset xtics
set xtics 0,30,179
set label 1 "(b)"
set tics scale 2
# set key default left bottom #font ",10"
set key at 110, graph 0.3 #font ",10"
# set key left bottom
plot "./data/ProtonPolDeg_100mev.dat" u 1:6  w line title "400 MeV" ls 6,\
	"" u 1:7  w line title "450 MeV"  ls 5,\
	"" u 1:8  w line title "500 MeV"  ls 3,\
	"" u 1:9  w line title "550 MeV"  ls 4,\
	"" u 1:11  w line title "AV18" ls 7

set label 1 "(c)"
set key at 110, graph 0.25 font ",8"
set xtics 0,30,180
# set autoscale

set style fill transparent solid 0.25 # partial transparency
# set style fill noborder # no separate top/bottom lines
plot "./dataSCS/Trunc_VU100MeV_Py_scs.dat" u 1:5 notitle lw 0.0001  lc rgb "#ff0000"   w filledcurves,\
	"./data/Trunc_VU100MeV_Py_sms.dat" u 1:6  notitle lw 0.0001 lc "blue"  w filledcurves,\
      "./dataSCS/Py_100mev_scs.dat" u 1:6 w l dt 2 lw 1 lc rgb "#DB0000" notitle,\
	"./data/ProtonPolDeg_100mev.dat" u 1:10  w l notitle ls 6 lw 1 lc "blue",\
	NaN w l title "SMS N^4LO+" ls 6 lw 1 lc "blue",\
	NaN title "SCS N^4LO" w l dt 2 lw 1 lc rgb "#DB0000"
	
