set term pdf enhanced size 5.00in, 2.70in
set output 'Cross_He_cutnum450.pdf'

MP_LEFT = .1
MP_RIGHT = .95
MP_BOTTOM = .11
MP_TOP = .99
MP_GAP = 0.04


set multiplot layout 2,4 margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP
#title "Chiral order dependence" \
              # tmargin 1
              # margins screen  0.1, 0.1, 0.5, 1 spacing screen MP_GAP
# set bmargin 5
# set lmargin 7
# set rmargin 0
# set tmargin 1
# set bmargin 0
  set style line 1 lt rgb "orange" lw 1.5 dt 1 pt 13 pi -25 ps 0.3 ## LO
  set style line 2 lt rgb "blue" lw 1.5 dt 1  pt 7 pi -25 ps 0.3 ## NLO
  set style line 3 lc rgb "green" lw 1.5 dt 2 ## N2LO
  set style line 4 lc rgb "red" lw 1.5 dt 3 ## N3LO
  set style line 5 lc rgb "black" lw 1.5 dt 1 ## N4LO
  set style line 6 lc rgb "cyan" lw 1.5 dt 5  ## N4LO+
  set style line 7 lc rgb "violet" lw 1.5  dt 4 ## AV18
  set pointintervalbox 0.02  ## interval to a point


set size square
# set yrange [0:40]
set ylabel "d^3{/Symbol s}/d{/Symbol W}_pdE_p [{/Symbol m}b sr^{-1} MeV^{-1}]" offset 2.5,-5 font ",11"
set xlabel "E_p [MeV]" font ",11"  offset 0,1,0
# set nokey
# unset key
set tics scale 2
set ytics font ",11" offset 0.6,0,0
set ytics 0.2
set xtics font ",11" offset 0,0.4,0
set xtics 30
set mytics 2
set mxtics 3
# set format x ""
# unset xtics

# set label 2 "d^2{/Symbol s}/d{/Symbol W}_pd{E_p} [{/Symbol m}b sr^{-1}]" at screen 0.03,0.9 #font ",12"


#################E = 120 MeV
################# CHIRAL

set yrange [0:0.8]
set key Left font ",10" at graph 0.8,graph 0.73 samplen 3

# set title "{{/Symbol q} = 0 deg}" font ",8"  offset 0,-0.5,0
set title "{{/Symbol q} = 0 deg}" font ",10"  offset 0,-2.5,0
plot "./data/ppn-CUTNUM450-120.0-MeV_LO-incl-0deg" every ::8 u 2:($3*10000)  w lp title "LO" ls 1,\
	"./data/ppn-CUTNUM450-120.0-MeV_NLO-incl-0deg" every ::8 u 2:($3*10000)  w lp title "NLO" ls 2,\
	"./data/ppn-CUTNUM450-120.0-MeV_N2LO-incl-0deg" every ::8 u 2:($3*10000)  w l notitle "N^2LO" ls 3,\
	"./data/ppn-CUTNUM450-120.0-MeV_N3LO-incl-0deg" every ::8 u 2:($3*10000)  w l notitle "N^3LO" ls 4,\
	"./data/ppn-CUTNUM450-120.0-MeV_N4LO-incl-0deg" every ::8 u 2:($3*10000)  w l notitle "N^4LO" ls 5,\
	"./data/ppn-CUTNUM450-120.0-MeV_N4LO+-incl-0deg" every ::8 u 2:($3*10000)  w l notitle "N^4LO+" ls 6

set autoscale	
unset ylabel
# unset ytics

# set format y ""

set title "{{/Symbol q} = 60 deg}"
plot "./data/ppn-CUTNUM450-120.0-MeV_LO-incl-60deg" every ::8 u 2:($3*10000)  w lp notitle "LO" ls 1,\
	"./data/ppn-CUTNUM450-120.0-MeV_NLO-incl-60deg" every ::8 u 2:($3*10000)  w lp notitle "NLO" ls 2,\
	"./data/ppn-CUTNUM450-120.0-MeV_N2LO-incl-60deg" every ::8 u 2:($3*10000)  w l title "N^2LO" ls 3,\
	"./data/ppn-CUTNUM450-120.0-MeV_N3LO-incl-60deg" every ::8 u 2:($3*10000)  w l title "N^3LO" ls 4,\
	"./data/ppn-CUTNUM450-120.0-MeV_N4LO-incl-60deg" every ::8 u 2:($3*10000)  w l notitle "N^4LO" ls 5,\
	"./data/ppn-CUTNUM450-120.0-MeV_N4LO+-incl-60deg" every ::8 u 2:($3*10000)  w l notitle "N^4LO+" ls 6

set ytics 0.2
set yrange [0:0.6]
set xrange [0:70]
set xtics 20
set key Left font ",10" at graph 0.93,graph 0.73

set title "{{/Symbol q} = 120 deg}"
plot "./data/ppn-CUTNUM450-120.0-MeV_LO-incl-120deg" every ::8 u 2:($3*10000)  w lp notitle "LO" ls 1,\
	"./data/ppn-CUTNUM450-120.0-MeV_NLO-incl-120deg" every ::8 u 2:($3*10000)  w lp notitle "NLO" ls 2,\
	"./data/ppn-CUTNUM450-120.0-MeV_N2LO-incl-120deg" every ::8 u 2:($3*10000)  w l notitle "N^2LO" ls 3,\
	"./data/ppn-CUTNUM450-120.0-MeV_N3LO-incl-120deg" every ::8 u 2:($3*10000)  w l notitle "N^3LO" ls 4,\
	"./data/ppn-CUTNUM450-120.0-MeV_N4LO-incl-120deg" every ::8 u 2:($3*10000)  w l title "N^4LO" ls 5,\
	"./data/ppn-CUTNUM450-120.0-MeV_N4LO+-incl-120deg" every ::8 u 2:($3*10000)  w l title "N^4LO+" ls 6

# set ytics 0,0.1,0.5
set ytics 0.2
set yrange [0:0.5]
set xrange [0:60]
unset key

set title "{{/Symbol q} = 180 deg}"
plot "./data/ppn-CUTNUM450-120.0-MeV_LO-incl-180deg" every ::8 u 2:($3*10000)  w lp title "LO" ls 1,\
	"./data/ppn-CUTNUM450-120.0-MeV_NLO-incl-180deg" every ::8 u 2:($3*10000)  w lp title "NLO" ls 2,\
	"./data/ppn-CUTNUM450-120.0-MeV_N2LO-incl-180deg" every ::8 u 2:($3*10000)  w l title "N^2LO" ls 3,\
	"./data/ppn-CUTNUM450-120.0-MeV_N3LO-incl-180deg" every ::8 u 2:($3*10000)  w l title "N^3LO" ls 4,\
	"./data/ppn-CUTNUM450-120.0-MeV_N4LO-incl-180deg" every ::8 u 2:($3*10000)  w l title "N^4LO" ls 5,\
	"./data/ppn-CUTNUM450-120.0-MeV_N4LO+-incl-180deg" every ::8 u 2:($3*10000)  w l title "N^4LO+" ls 6

set autoscale
# set ytics autofreq
set ytics 0.2

#################E = 120 MeV
################# CUT-OFF

# set key font ",6" top right
set key Left font ",10" at graph 1.0,graph 0.73
set xtics 30
set title "{{/Symbol q} = 0 deg}"
plot "./data/ppn-CUTNUM400-120.0-MeV_N4LO-incl-0deg" every ::9 u 2:($3*10000)  w l ls 6 title "400 MeV",\
	"./data/ppn-CUTNUM450-120.0-MeV_N4LO-incl-0deg" every ::9 u 2:($3*10000)  w line ls 5 title "450 MeV", \
	"./data/ppn-CUTNUM500-120.0-MeV_N4LO-incl-0deg" every ::9 u 2:($3*10000)  w line ls 3 notitle "500 MeV",\
	"./data/ppn-CUTNUM550-120.0-MeV_N4LO-incl-0deg" every ::9 u 2:($3*10000)  w line ls 4 notitle "550 MeV"

# unset ylabel
# unset ytics
# set format y ""

# set title "{{/Symbol q} = 60 deg}" font ",8"  offset 0,-0.5,0
set title "{{/Symbol q} = 60 deg}"
plot "./data/ppn-CUTNUM400-120.0-MeV_N4LO-incl-60deg" every ::9 u 2:($3*10000)  w l ls 6 notitle "{400 MeV",\
	"./data/ppn-CUTNUM450-120.0-MeV_N4LO-incl-60deg" every ::9 u 2:($3*10000)  w line ls 5 notitle "{450 MeV", \
	"./data/ppn-CUTNUM500-120.0-MeV_N4LO-incl-60deg" every ::9 u 2:($3*10000)  w line ls 3 title "{500 MeV",\
	"./data/ppn-CUTNUM550-120.0-MeV_N4LO-incl-60deg" every ::9 u 2:($3*10000)  w line ls 4 title "{550 MeV"

unset key


set xrange [0:70]
set xtics 20
set title "{{/Symbol q} = 120 deg}"
plot "./data/ppn-CUTNUM400-120.0-MeV_N4LO-incl-120deg" every ::9 u 2:($3*10000)  w l ls 6 title "{/Symbol L} = 400 MeV",\
	"./data/ppn-CUTNUM450-120.0-MeV_N4LO-incl-120deg" every ::9 u 2:($3*10000)  w line ls 5 title "{/Symbol L} = 450 MeV", \
	"./data/ppn-CUTNUM500-120.0-MeV_N4LO-incl-120deg" every ::9 u 2:($3*10000)  w line ls 3 title "{/Symbol L} = 500 MeV",\
	"./data/ppn-CUTNUM550-120.0-MeV_N4LO-incl-120deg" every ::9 u 2:($3*10000)  w line ls 4 title "{/Symbol L} = 500 MeV"

# set ytics 0,0.1,0.5
set xrange [0:60]
set yrange [0:0.5]
set ytics 0.2
set title "{{/Symbol q} = 180 deg}"
plot "./data/ppn-CUTNUM400-120.0-MeV_N4LO-incl-180deg" every ::9 u 2:($3*10000)  w l ls 6 title "{/Symbol L} = 400 MeV",\
	"./data/ppn-CUTNUM450-120.0-MeV_N4LO-incl-180deg" every ::9 u 2:($3*10000)  w line ls 5 title "{/Symbol L} = 450 MeV", \
	"./data/ppn-CUTNUM500-120.0-MeV_N4LO-incl-180deg" every ::9 u 2:($3*10000)  w line ls 3 title "{/Symbol L} = 500 MeV",\
	"./data/ppn-CUTNUM550-120.0-MeV_N4LO-incl-180deg" every ::9 u 2:($3*10000)  w line ls 4 title "{/Symbol L} = 500 MeV"

