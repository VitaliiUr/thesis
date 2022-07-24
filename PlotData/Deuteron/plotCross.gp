set term pdf enhanced
set output 'Cross_All.pdf'

MP_LEFT = .1
MP_RIGHT = .95
MP_BOTTOM = .14
MP_TOP = .95
MP_GAP = 0.0


set multiplot layout 2,3 \
				margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP
              # margins screen  0.1, 0.1, 0.5, 1 spacing screen MP_GAP
              # tmargin 1

# set lmargin 7
# set rmargin 0
# set tmargin 1
# set bmargin 0
  set style line 1 lt rgb "orange" lw 1.5 dt 1 pt 13 pi -10 ps 0.3 ## LO
  set style line 2 lt rgb "blue" lw 1.5 dt 1  pt 7 pi -10 ps 0.3 ## NLO
  set style line 3 lc rgb "green" lw 1.5 dt 2 ## N2LO
  set style line 4 lc rgb "red" lw 1.5 dt 3 ## N3LO
  set style line 5 lc rgb "black" lw 1.5 dt 1 ## N4LO
  set style line 6 lc rgb "cyan" lw 1.5 dt 5  ## N4LO+
  set style line 7 lc rgb "violet" lw 1.5  dt 4 ## AV18
  set pointintervalbox 0.02  ## interval to a point

# set size square
set tics scale 2
set yrange [3:40]
set ylabel "d^2{/Symbol s}/d{/Symbol W} [{/Symbol m}b sr^{-1}]" offset 1.2, -4.7
set key at 117,25 font ",8"
# unset key
set ytics 0,8,40
set xtics 0,30,170
set mytics 4
set mxtics 3
set format x ""
# unset xtics

set label 1 "(a)" at 130, graph 0.8

# set label "d^2{/Symbol s}/d{/Symbol W} [{/Symbol m}b sr^{-1}]" at screen 0.03,0.4 rotate by 90

plot "./data/CROSS_30mev.dat" u 1:2  w lp title "LO" ls 1,\
	"./data/CROSS_30mev.dat" u 1:3  w lp title "NLO" ls 2,\
	"./data/CROSS_30mev.dat" u 1:4  w l title "N^2LO" ls 3,\
	"./data/CROSS_30mev.dat" u 1:5  w l title "N^3LO" ls 4,\
	"./data/CROSS_30mev.dat" u 1:7  w l title "N^4LO" ls 5,\
	"./data/CROSS_30mev.dat" u 1:10  w l title "N^4LO+" ls 6,\
	 "./data/CROSS_30mev.dat" u 1:11  w l title "AV18" ls 7,\
	"./ExpData/Exp30.dat" u 1:2:3 notitle pt 7 ps 0.5 lc rgb "black" with yerrorbars,\
	 "./ExpData/Exp30_2.dat" u 1:2:3 notitle pt 6 ps 0.5 lc rgb "black" with yerrorbars

unset ylabel
# unset ytics
set format y ""
set key at 120,22 font ",8"
set label 1 "(b)"

plot "./data/Trunc_VU30MeV_cross_sms.dat" u 1:2 lc rgb "#ffeda0" lw 0.0001  w filledcurves title "NLO",\
	"./data/Trunc_VU30MeV_cross_sms.dat" u 1:3 lc rgb "#a1d76a" lw 0.0001  w filledcurves title "N^2LO", \
	"./data/Trunc_VU30MeV_cross_sms.dat" u 1:4  lc rgb "#3399ff" lw 0.0001  w filledcurves title "N^3LO",\
	"./data/Trunc_VU30MeV_cross_sms.dat" u 1:5  title "N^4LO" lw 0.0001  lc rgb "#ff0000" w filledcurves, \
	"./data/Trunc_VU30MeV_cross_sms.dat" u 1:6  title "N^4LO+" lw 0.0001  lc rgb "black" w filledcurves,\
	"./ExpData/Exp30.dat" u 1:2:3 notitle pt 7 ps 0.5 lc rgb "black" with yerrorbars,\
	 "./ExpData/Exp30_2.dat" u 1:2:3 notitle pt 6 ps 0.5 lc rgb "black" with yerrorbars

set key at 124,19 font ",8"
set label 1 "(c)"

plot "./data/CROSS_30mev.dat" u 1:6  w line title "400 MeV" ls 6,\
	"./data/CROSS_30mev.dat" u 1:7  w line title "450 MeV" ls 5,\
	"./data/CROSS_30mev.dat" u 1:8  w line title "500 MeV" ls 3,\
	"./data/CROSS_30mev.dat" u 1:9  w line title "550 MeV"  ls 4,\
	"./data/CROSS_30mev.dat" u 1:11  w line title "AV18" ls 7, \
	"./ExpData/Exp30.dat" u 1:2:3 notitle pt 7 ps 0.5 lc rgb "black" with yerrorbars,\
	"./ExpData/Exp30_2.dat" u 1:2:3 notitle pt 6 ps 0.5 lc rgb "black" with yerrorbars

# set tics scale 2
set yrange [1:8]
set format x "%1.f"
set format y "%1.f"
set ytics 1,2,7.9
# set xtics 0,30,170
set mytics 4
set mxtics 3
unset key
set xlabel " {/Symbol Q}_p [deg]" offset 0,0.4
set label 1 "(d)"

plot "./data/CROSS_100mev.dat" u 1:2  w lp title "LO" ls 1,\
"./data/CROSS_100mev.dat" u 1:3  w lp title "NLO" ls 2,\
	"./data/CROSS_100mev.dat" u 1:4  w l title "N^2LO" ls 3,\
	"./data/CROSS_100mev.dat" u 1:5  w l title "N^3LO" ls 4,\
	"./data/CROSS_100mev.dat" u 1:7  w l title "N^4LO" ls 5,\
	"./data/CROSS_100mev.dat" u 1:10  w l title "N^4LO+" ls 6,\
	"./data/CROSS_100mev.dat" u 1:11  w l title "AV18" ls 7,\
	"./ExpData/Exp100.dat" u 1:2:3 notitle pt 7 ps 0.5 lc rgb "black" with yerrorbars,\
	"./ExpData/Exp100_2.dat" u 1:2:3 notitle pt 6 ps 0.5 lc rgb "black" with yerrorbars,\
	"./ExpData/Exp100n.dat" u 1:2:3 notitle pt 5 ps 0.5 lc rgb "black" with yerrorbars,\
	"./ExpData/Exp100_3.dat" u 1:2:3 notitle pt 9 ps 0.5 lc rgb "black" with yerrorbars

set format y ""
set label 1 "(e)"

plot "./data/Trunc_VU100MeV_cross_sms.dat" u 1:2  w filledcurves lc rgb "#ffeda0" lw 0.0001 title "NLO",\
	"./data/Trunc_VU100MeV_cross_sms.dat" u 1:3  lc rgb "#a1d76a" lw 0.0001  w filledcurves title "N^2LO",\
	"./data/Trunc_VU100MeV_cross_sms.dat" u 1:4  lc rgb "#3399ff" lw 0.0001  w filledcurves title "N^3LO",\
	"./data/Trunc_VU100MeV_cross_sms.dat" u 1:5  lc rgb "#ff0000" lw 0.0001  title "N^4LO" w filledcurves, \
	"./data/Trunc_VU100MeV_cross_sms.dat" u 1:6  lc rgb "black" lw 0.0001  title "N^4LO+" w filledcurves, \
	"./ExpData/Exp100.dat" u 1:2:3 notitle pt 7 ps 0.5 lc rgb "black" with yerrorbars,\
	"./ExpData/Exp100_2.dat" u 1:2:3 notitle pt 6 ps 0.5 lc rgb "black" with yerrorbars,\
	"./ExpData/Exp100n.dat" u 1:2:3 notitle pt 5 ps 0.5 lc rgb "black" with yerrorbars,\
	"./ExpData/Exp100_3.dat" u 1:2:3 notitle pt 9 ps 0.5 lc rgb "black" with yerrorbars

set xtics 0,30,180
set label 1 "(f)"

plot "./data/CROSS_100mev.dat" u 1:6  w line title "400 MeV" ls 6,\
	"./data/CROSS_100mev.dat" u 1:7  w line title "450 MeV" ls 5,\
	"./data/CROSS_100mev.dat" u 1:8  w line title "500 MeV" ls 3,\
	"./data/CROSS_100mev.dat" u 1:9  w line title "550 MeV"  ls 4,\
	"./data/CROSS_100mev.dat" u 1:11  w line title "AV18" ls 7, \
	"./ExpData/Exp100.dat" u 1:2:3 notitle pt 7 ps 0.5 lc rgb "black" with yerrorbars,\
	"./ExpData/Exp100_2.dat" u 1:2:3 notitle pt 6 ps 0.5 lc rgb "black" with yerrorbars,\
	"./ExpData/Exp100n.dat" u 1:2:3 notitle pt 5 ps 0.5 lc rgb "black" with yerrorbars,\
	"./ExpData/Exp100_3.dat" u 1:2:3 notitle pt 9 ps 0.5 lc rgb "black" with yerrorbars

