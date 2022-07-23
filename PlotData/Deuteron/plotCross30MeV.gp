set term pdf enhanced
set output 'Cross_All.pdf'

MP_LEFT = .1
MP_RIGHT = .95
MP_BOTTOM = .14
MP_TOP = .9
MP_GAP = 0.0


set multiplot layout 2,3 \
				margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP
              # margins screen  0.1, 0.1, 0.5, 1 spacing screen MP_GAP
              # tmargin 1

# set lmargin 7
# set rmargin 0
# set tmargin 1
# set bmargin 0
  set style line 1 lt rgb "orange" lw 2.5 dt 1 #pt 6 pi -10 ps 0.5
  set style line 2 lt rgb "blue" lw 2.5 dt 4
  set style line 3 lc rgb "green" lw 2.5 dt 2
  set style line 4 lc rgb "red" lw 2.5 dt 3
  set style line 5 lc rgb "black" lw 2.5 dt 1
  set style line 6 lc rgb "violet" lw 2.5 dt 5
  set style line 7 lc rgb "cyan" lw 2.5
  set pointintervalbox 0.02  ## interval to a point

# set size square
set yrange [0:40]
# set ylabel "d^2{/Symbol s}/d{/Symbol W} [{/Symbol m}b sr^{-1}]"
set key center bottom font ",8"
set ytics 0,8,40
set mytics 2
set mxtics 2
set format x ""
unset xtics

set label 2 "d^2{/Symbol s}/d{/Symbol W} [{/Symbol m}b sr^{-1}]" at screen 0.03,0.4 rotate by 90

plot "CROSS_30mev.dat" u 1:2  w l title "LO" ls 1,"CROSS_30mev.dat" u 1:3  w line title "NLO" ls 2,"CROSS_30mev.dat" u 1:4  w line title "N^2LO" ls 3,"CROSS_30mev.dat" u 1:5  w line title "N^3LO" ls 4,"CROSS_30mev.dat" u 1:7  w line title "N^4LO" ls 5,"CROSS_30mev.dat" u 1:10  w line title "N^4LO+" ls 7, "CROSS_30mev.dat" u 1:11  w line title "AV18" ls 6, \
	"Exp30.dat" u 1:2:3 notitle pt 7 ps 0.5 lc rgb "black" with yerrorbars, "Exp30_2.dat" u 1:2:3 notitle pt 6 ps 0.5 lc rgb "black" with yerrorbars

unset ylabel
unset ytics
set format y ""

plot "Thruns_VU30MeV.dat" u 1:2 lc rgb "#ffeda0"  w filledcurves title "NLO","Thruns_VU30MeV.dat" u 1:3 lc rgb "#a1d76a"  w filledcurves title "N^2LO", \
	"Thruns_VU30MeV.dat" u 1:4  lc rgb "#3399ff"  w filledcurves title "N^3LO","Thruns_VU30MeV.dat" u 1:5  title "N^4LO"  lc rgb "#ff0000" w filledcurves, \
	"Exp30.dat" u 1:2:3 notitle pt 7 ps 0.5 lc rgb "black" with yerrorbars, "Exp30_2.dat" u 1:2:3 notitle pt 6 ps 0.5 lc rgb "black" with yerrorbars


plot "CROSS_30mev.dat" u 1:6  w line title "{/Symbol L} = 400 MeV" ls 1, "CROSS_30mev.dat" u 1:7  w line title "{/Symbol L} = 450 MeV" ls 2, "CROSS_30mev.dat" u 1:8  w line title "{/Symbol L} = 500 MeV" ls 3, "CROSS_30mev.dat" u 1:9  w line title "{/Symbol L} = 550 MeV"  ls 4, "CROSS_30mev.dat" u 1:11  w line title "AV18" ls 6, \
	"Exp30.dat" u 1:2:3 notitle pt 7 ps 0.5 lc rgb "black" with yerrorbars, "Exp30_2.dat" u 1:2:3 notitle pt 6 ps 0.5 lc rgb "black" with yerrorbars


set yrange [0:8]
set mytics 2
set format x "%1.f"
set format y "%1.f"
set ytics 0,2,7
set xtics 0,30,170
unset key
set xlabel " {/Symbol Q}_p [deg]" offset 0,0.4

plot "CROSS_100mev.dat" u 1:2  w l title "LO" ls 1,"CROSS_100mev.dat" u 1:3  w line title "NLO" ls 2, \
	"CROSS_100mev.dat" u 1:4  w line title "N^2LO" ls 3,"CROSS_100mev.dat" u 1:5  w line title "N^3LO" ls 4,"CROSS_100mev.dat" u 1:7  w line title "N^4LO" ls 5,"CROSS_100mev.dat" u 1:10  w line title "N^4LO+" ls 7, "CROSS_100mev.dat" u 1:11  w line title "AV18" ls 6, \
	"Exp100.dat" u 1:2:3 notitle pt 7 ps 0.5 lc rgb "black" with yerrorbars, "Exp100_2.dat" u 1:2:3 notitle pt 6 ps 0.5 lc rgb "black" with yerrorbars, "Exp100n.dat" u 1:2:3 notitle pt 5 ps 0.5 lc rgb "black" with yerrorbars

set format y ""

plot "Thruns_VU100MeV.dat" u 1:2  w filledcurves lc rgb "#ffeda0" title "NLO","Thruns_VU100MeV.dat" u 1:3  lc rgb "#a1d76a"  w filledcurves title "N^2LO", \
	"Thruns_VU100MeV.dat" u 1:4  lc rgb "#3399ff"  w filledcurves title "N^3LO","Thruns_VU100MeV.dat" u 1:5  lc rgb "#ff0000"  title "N^4LO" w filledcurves, \
	"Exp100.dat" u 1:2:3 notitle pt 7 ps 0.5 lc rgb "black" with yerrorbars, "Exp100_2.dat" u 1:2:3 notitle pt 6 ps 0.5 lc rgb "black" with yerrorbars, "Exp100n.dat" u 1:2:3 notitle pt 5 ps 0.5 lc rgb "black" with yerrorbars

set xtics 0,30,180

plot "CROSS_100mev.dat" u 1:6  w line title "{/Symbol L} = 400 MeV" ls 1, "CROSS_100mev.dat" u 1:7  w line title "{/Symbol L} = 450 MeV" ls 2, "CROSS_100mev.dat" u 1:8  w line title "{/Symbol L} = 500 MeV" ls 3, "CROSS_100mev.dat" u 1:9  w line title "{/Symbol L} = 550 MeV"  ls 4, "CROSS_100mev.dat" u 1:11  w line title "AV18" ls 6, \
	"Exp100.dat" u 1:2:3 notitle pt 7 ps 0.5 lc rgb "black" with yerrorbars, "Exp100_2.dat" u 1:2:3 notitle pt 6 ps 0.5 lc rgb "black" with yerrorbars, "Exp100n.dat" u 1:2:3 notitle pt 5 ps 0.5 lc rgb "black" with yerrorbars

