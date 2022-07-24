set term pdfcairo enhanced size 5.00in, 2.70in

MP_LEFT = .1
MP_RIGHT = .95
MP_BOTTOM = .17
MP_TOP = .95
MP_GAP = 0.0005

  set style line 1 lt rgb "orange" lw 2 dt 1 pt 13 pi -10 ps 0.3 ## LO
  set style line 2 lt rgb "blue" lw 2 dt 1  pt 7 pi -10 ps 0.3 ## NLO
  set style line 3 lc rgb "green" lw 2 dt 2 ## N2LO
  set style line 4 lc rgb "red" lw 2 dt 3 ## N3LO
  set style line 5 lc rgb "black" lw 2 dt 1 ## N4LO
  set style line 6 lc rgb "cyan" lw 2 dt 5  ## N4LO+
  set style line 7 lc rgb "violet" lw 2  dt 4 ## AV18
  set pointintervalbox 0.02  ## interval to a point


set output 'AX_100MeV_all.pdf'
set multiplot layout 1,2 \
              margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP

set tics scale 2
set key at 100, graph 0.42 font ",10"
# unset key
set ylabel 'A_X' offset 0.9,0
set xlabel '{/Symbol q}_p [deg]' 
# set format y "%g"
set xtics 0,30,179
set mxtics 3
set mytics 2
set yrange [0:0.65]
set label 1 "(a)" at 130, graph 0.8 font ",14"

plot "./data/DeutAX_Deg_100mev.dat" u 1:(-$2)  w lp title "LO" ls 1,\
	"./data/DeutAX_Deg_100mev.dat" u 1:(-$3)  w lp title "NLO" ls 2,\
	"./data/DeutAX_Deg_100mev.dat" u 1:(-$4)  w l title "N^2LO" ls 3,\
	"./data/DeutAX_Deg_100mev.dat" u 1:(-$5)  w l title "N^3LO" ls 4,\
	"./data/DeutAX_Deg_100mev.dat" u 1:(-$7)  w l title "N^4LO" ls 5,\
	"./data/DeutAX_Deg_100mev.dat" u 1:(-$10)  w l title "N^4LO+" ls 6,\
	"./data/DeutAX_Deg_100mev.dat" u 1:(-$11)  w l title "AV18" ls 7

unset ylabel
set format y ""
unset xtics
set tics scale 2
# set key bottom center
# set key font ",12"
# unset key
# set ytics nomirror
unset ytics

set key at 110, graph 0.32 font ",10"

set y2tics autofreq
set y2range [0:0.65]
set y2tics mirror
set my2tics 2
set y2tics 0.1
set xtics 0,30,180
set label 1 "(b)"
set format y2 ""
#set key bottom center font ",16"
plot "./data/DeutAX_Deg_100mev.dat" u 1:(-$6)  w line title "400 MeV" ls 6  axis x1y2,\
	"./data/DeutAX_Deg_100mev.dat" u 1:(-$7)  w line title "450 MeV" ls 5  axis x1y2,\
	"./data/DeutAX_Deg_100mev.dat" u 1:(-$8)  w line title "500 MeV" ls 3  axis x1y2,\
	"./data/DeutAX_Deg_100mev.dat" u 1:(-$9)  w line title "550 MeV"  ls 4  axis x1y2,\
	"./data/DeutAX_Deg_100mev.dat" u 1:(-$11)  w line title "AV18" ls 7 axis x1y2
unset multiplot

# set output 'AX_30MeV_all.pdf'
# set multiplot layout 1,2 \
#               margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP


# set key top center
# set ylabel 'A_X' offset 0.9,0
# set xlabel '{/Symbol q}_p [deg]' 
# set format y "%g"
# set xtics 0,20,170
# set ytics 0.1
# set mxtics 4
# set mytics 4
# set yrange [-0.9:0]
# set label 1 "(a)" at 20, graph 0.12

# plot "./data/DeutAX_Deg_30mev.dat" u 1:2  w lp title "LO" ls 1,\
# 	"./data/DeutAX_Deg_30mev.dat" u 1:3  w lp title "NLO" ls 2,\
# 	"./data/DeutAX_Deg_30mev.dat" u 1:4  w l title "N^2LO" ls 3,\
# 	"./data/DeutAX_Deg_30mev.dat" u 1:5  w l title "N^3LO" ls 4,\
# 	"./data/DeutAX_Deg_30mev.dat" u 1:7  w l title "N^4LO" ls 5,\
# 	"./data/DeutAX_Deg_30mev.dat" u 1:10  w l title "N^4LO+" ls 6,\
# 	"./data/DeutAX_Deg_30mev.dat" u 1:11  w l title "AV18" ls 7

# unset ylabel
# set format y ""
# unset xtics
# set xtics 0,20,180
# #unset ytics

# set key font ",12"
# set label 1 "(b)"
# #set key bottom center font ",16"
# plot "./data/DeutAX_Deg_30mev.dat" u 1:6  w line title "400 MeV" ls 6,\
# 	"./data/DeutAX_Deg_30mev.dat" u 1:7  w line title "450 MeV" ls 5,\
# 	"./data/DeutAX_Deg_30mev.dat" u 1:8  w line title "500 MeV" ls 3,\
# 	"./data/DeutAX_Deg_30mev.dat" u 1:9  w line title "550 MeV"  ls 4,\
# 	"./data/DeutAX_Deg_30mev.dat" u 1:11  w line title "AV18" ls 7
# unset multiplot
