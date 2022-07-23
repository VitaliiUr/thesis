set term pdf enhanced size 5.00in, 5.00in
set output 'Cross_chiral_orders_120MeV.pdf'

MP_LEFT = .1
MP_RIGHT = .95
MP_BOTTOM = .07
MP_TOP = .92
MP_GAP = 0.0


set multiplot layout 2,2 title "Chiral order dependence" \
				margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP
              # tmargin 1
              # margins screen  0.1, 0.1, 0.5, 1 spacing screen MP_GAP
# set bmargin 5
# set lmargin 7
# set rmargin 0
# set tmargin 1
# set bmargin 0
  set style line 1 lt rgb "orange" lw 1.5 dt 1 #pt 6 pi -10 ps 0.5
  set style line 2 lt rgb "blue" lw 1.5 dt 4
  set style line 3 lc rgb "green" lw 1.5 dt 2
  set style line 4 lc rgb "red" lw 1.5 dt 3
  set style line 5 lc rgb "black" lw 1.5 dt 1
  set style line 6 lc rgb "violet" lw 1.5 dt 5
  set style line 7 lc rgb "cyan" lw 1.5
  set pointintervalbox 0.02  ## interval to a point

set size square
# set yrange [0:40]
set ylabel "d^2{/Symbol s}/d{/Symbol W} [{/Symbol m}b sr^{-1}]" offset 2,-8 font ",10"
# set xlabel "E_p [MeV]" font ",8"  offset 0,1,0
# set nokey
set key font ",10"
set yrange [0:0.8]
set ytics 0,0.1,0.8 font ",6"
set xtics font ",6" offset 0,1.7,0
set mytics 2
# set mxtics 2
# set format x ""
# unset xtics

# set label "d^2{/Symbol s}/d{/Symbol W}_pd{E_p} [{/Symbol m}b sr^{-1}]" at screen 0.03,0.9 font ",10" rotate by 90

#################E = 120 MeV

set title "{{/Symbol q} = 0 deg}" font ",10"  offset 0,-2.1,0
plot "./data/ppn-120.0-MeV_LO-incl-0deg" every ::8 u 2:($3*10000)  w l title "LO" ls 1,\
	"./data/ppn-120.0-MeV_NLO-incl-0deg" every ::8 u 2:($3*10000)  w line title "NLO" ls 2,\
	"./data/ppn-120.0-MeV_N2LO-incl-0deg" every ::8 u 2:($3*10000)  w line title "N^2LO" ls 3,\
	"./data/ppn-120.0-MeV_N3LO-incl-0deg" every ::8 u 2:($3*10000)  w line title "N^3LO" ls 4,\
	"./data/ppn-120.0-MeV_N4LO-incl-0deg" every ::8 u 2:($3*10000)  w line title "N^4LO" ls 5,\
	"./data/ppn-120.0-MeV_N4LO+-incl-0deg" every ::8 u 2:($3*10000)  w line title "N^4LO+" ls 6
unset ylabel
# unset ytics
# set y2tics font ",6"
unset key
set format y ""

set title "{{/Symbol q} = 60 deg}"
plot "./data/ppn-120.0-MeV_LO-incl-60deg" every ::8 u 2:($3*10000)  w l title "LO" ls 1,\
	"./data/ppn-120.0-MeV_NLO-incl-60deg" every ::8 u 2:($3*10000)  w line title "NLO" ls 2,\
	"./data/ppn-120.0-MeV_N2LO-incl-60deg" every ::8 u 2:($3*10000)  w line title "N^2LO" ls 3,\
	"./data/ppn-120.0-MeV_N3LO-incl-60deg" every ::8 u 2:($3*10000)  w line title "N^3LO" ls 4,\
	"./data/ppn-120.0-MeV_N4LO-incl-60deg" every ::8 u 2:($3*10000)  w line title "N^4LO" ls 5,\
	"./data/ppn-120.0-MeV_N4LO+-incl-60deg" every ::8 u 2:($3*10000)  w line title "N^4LO+" ls 6

# unset y2tics
set yrange [0:0.6]
set ytics 0,0.05,0.58 font ",6"
set xlabel "E_p [MeV]" font ",8"  offset 0,1,0
set xtics font ",6" offset 0,0.7,0
set format x "%g"
set format y "%g"
# set ylabel "d^2{/Symbol s}/d{/Symbol W} [{/Symbol m}b sr^{-1}]" offset 2,0 font ",10"

# set ytics font ",6" offset 0.5,0,0

set title "{{/Symbol q} = 120 deg}"
plot "./data/ppn-120.0-MeV_LO-incl-120deg" every ::8 u 2:($3*10000)  w l title "LO" ls 1,\
	"./data/ppn-120.0-MeV_NLO-incl-120deg" every ::8 u 2:($3*10000)  w line title "NLO" ls 2,\
	"./data/ppn-120.0-MeV_N2LO-incl-120deg" every ::8 u 2:($3*10000)  w line title "N^2LO" ls 3,\
	"./data/ppn-120.0-MeV_N3LO-incl-120deg" every ::8 u 2:($3*10000)  w line title "N^3LO" ls 4,\
	"./data/ppn-120.0-MeV_N4LO-incl-120deg" every ::8 u 2:($3*10000)  w line title "N^4LO" ls 5,\
	"./data/ppn-120.0-MeV_N4LO+-incl-120deg" every ::8 u 2:($3*10000)  w line title "N^4LO+" ls 6


# unset ytics
# set y2tics 0,0.05, 0.58 font ",6"
set format y ""
unset ylabel


set title "{{/Symbol q} = 180 deg}"
plot "./data/ppn-120.0-MeV_LO-incl-180deg" every ::8 u 2:($3*10000)  w l title "LO" ls 1,\
	"./data/ppn-120.0-MeV_NLO-incl-180deg" every ::8 u 2:($3*10000)  w line title "NLO" ls 2,\
	"./data/ppn-120.0-MeV_N2LO-incl-180deg" every ::8 u 2:($3*10000)  w line title "N^2LO" ls 3,\
	"./data/ppn-120.0-MeV_N3LO-incl-180deg" every ::8 u 2:($3*10000)  w line title "N^3LO" ls 4,\
	"./data/ppn-120.0-MeV_N4LO-incl-180deg" every ::8 u 2:($3*10000)  w line title "N^4LO" ls 5,\
	"./data/ppn-120.0-MeV_N4LO+-incl-180deg" every ::8 u 2:($3*10000)  w line title "N^4LO+" ls 6
