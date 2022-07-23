# set term pdf enhanced
# set term postscript eps enhanced color size 3.5,3.5/2
set term pdfcairo enhanced size 5.00in, 2.00in
set output 'Cross140MeV.pdf'

MP_LEFT = .1
MP_RIGHT = .95
MP_BOTTOM = .2
MP_TOP = 1.
MP_GAP = 0.0


set multiplot layout 1,3 \
		margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP\
				spacing screen MP_GAP
				# title "Cross-section 140 MeV"
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

set size square
set yrange [1:6]
set format x "%1.f"
set format y "%1.f"

set ylabel "d^2{/Symbol s}/d{/Symbol W} [{/Symbol m}b sr^{-1}]"
# set key center bottom font ",8"
set key right top font ",8"
set ytics 1,1,6
set mytics 2
set mxtics 2
# set format x ""
set xtics 0, 30, 170
set xlabel " {/Symbol Q}_p [deg]" offset 0,0.4

# set label 2 "d^2{/Symbol s}/d{/Symbol W} [{/Symbol m}b sr^{-1}]" at screen 0.03,0.4 rotate by 90

plot "./data/CROSS_140mev.dat" u 1:2  w lp title "LO" ls 1,\
	"./data/CROSS_140mev.dat" u 1:3  w lp title "NLO" ls 2,\
	"./data/CROSS_140mev.dat" u 1:4  w line title "N^2LO" ls 3,\
	"./data/CROSS_140mev.dat" u 1:5  w line title "N^3LO" ls 4,\
	"./data/CROSS_140mev.dat" u 1:7  w line title "N^4LO" ls 5,\
	"./data/CROSS_140mev.dat" u 1:10  w l title "N^4LO+" ls 7,\
	"./ExpData/Exp140.dat" u 1:2:3 notitle pt 7 ps 0.5 lc rgb "black" with yerrorbars
	# "CROSS_140mev.dat" u 1:11  w line title "AV18" ls 6, \

unset ylabel
unset ytics
set format y ""

plot "./data/Thruns_VU140MeV.dat" u 1:2 lc rgb "#ffeda0"  w filledcurves title "NLO",\
	"./data/Thruns_VU140MeV.dat" u 1:3 lc rgb "#a1d76a"  w filledcurves title "N^2LO", \
	"./data/Thruns_VU140MeV.dat" u 1:4  lc rgb "#3399ff"  w filledcurves title "N^3LO",\
	"./data/Thruns_VU140MeV.dat" u 1:5  title "N^4LO"  lc rgb "#ff0000" w filledcurves, \
	"./data/Thruns_VU140MeV.dat" u 1:6  title "N^4LO+"  lc rgb "black" w filledcurves,\
	"./ExpData/Exp140.dat" u 1:2:3 notitle pt 7 ps 0.5 lc rgb "black" with yerrorbars

set xtics 0, 30, 180

set key left bottom font ",8"

plot "./data/CROSS_140mev.dat" u 1:6  w line title "{/Symbol L} = 400 MeV" ls 6,\
	"./data/CROSS_140mev.dat" u 1:7  w line title "{/Symbol L} = 450 MeV" ls 5,\
	"./data/CROSS_140mev.dat" u 1:8  w line title "{/Symbol L} = 500 MeV" ls 3,\
	"./data/CROSS_140mev.dat" u 1:9  w line title "{/Symbol L} = 550 MeV"  ls 4,\
	"./ExpData/Exp140.dat" u 1:2:3 notitle pt 7 ps 0.5 lc rgb "black" with yerrorbars
 # "CROSS_140mev.dat" u 1:11  w line title "AV18" ls 6, \



unset multiplot
# set term postscript eps enhanced color size 2*0.9, 2*0.9
# set term eps enhanced color size 5*0.9, 5*0.9
# set term pdf enhanced size 2*0.9, 2*0.9
set term pdfcairo enhanced size 5.00in, 3.00in	
set output 'Cross140MeV_SingCur.pdf'
set format y "%1.f"
set yrange [0:6]
set ytics 0,1,6
set ylabel "d^2{/Symbol s}/d{/Symbol W} [{/Symbol m}b sr^{-1}]"
# set key right top font ",8"
unset key

set xtics 0, 30, 180

# set title "Cross-section 140 MeV N^4LO"
plot "./data/CROSS_SingCur_140mev.dat" u 1:2  w line title "Single-nucleon current" ls 1 lw 2.5,\
	"./data/CROSS_SingCur_140mev.dat" u 1:3  w line title "Siegert" ls 2 lw 2.5,\
	"./ExpData/Exp140.dat" u 1:2:3 notitle pt 7 ps 0.5 lc rgb "black" with yerrorbars