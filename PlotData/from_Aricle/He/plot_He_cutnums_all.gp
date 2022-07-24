set term pdf enhanced

MP_LEFT = .1
MP_RIGHT = .95
MP_BOTTOM = .1
MP_TOP = .85
MP_GAP = 0.05

set output 'Cross_N4LO_cutnums.pdf'
set multiplot layout 2,4 title "Cutnum dependence N4LO"\
				margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP
              # tmargin 1
              # margins screen  0.1, 0.1, 0.5, 1 spacing screen MP_GAP
# set bmargin 5
# set lmargin 7
# set rmargin 0
# set tmargin 1
# set bmargin 0
  set style line 1 lt rgb "orange" lw 1 dt 1 #pt 6 pi -10 ps 0.5
  set style line 2 lt rgb "blue" lw 1 dt 4
  set style line 3 lc rgb "green" lw 1 dt 2
  set style line 4 lc rgb "red" lw 1 dt 3
  set style line 5 lc rgb "black" lw 1 dt 1
  set style line 6 lc rgb "violet" lw 1 dt 5
  set style line 7 lc rgb "cyan" lw 1
  set pointintervalbox 0.02  ## interval to a point

set size square
# set yrange [0:40]
# set ylabel "d^2{/Symbol s}/d{/Symbol W} [{/Symbol m}b sr^{-1}]" font ",6"
set xlabel "E_p [MeV]" font ",6"  offset 0,1.6,0
set nokey
set ytics font ",6" offset 0.5,0,0
set xtics font ",6" offset 0,0.7,0
set mytics 2
set mxtics 2
# set format x ""
# unset xtics

set label 2 "d^2{/Symbol s}/d{/Symbol W}_pd{E_p} [{/Symbol m}b sr^{-1}]" at screen 0.03,0.9 font ",12"

#################E = 40 MeV

set title "{{/Symbol q} = 0 deg}" font ",6"  offset 0,-0.5,0
plot "./data/ppn-CUTNUM400-40.0-MeV_N4LO-incl-0deg-2020-01-22" every ::9 u 2:($3*10000)  w l ls 1 title "{/Symbol L} = 400 MeV",\
	"./data/ppn-CUTNUM450-40.0-MeV_N4LO-incl-0deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 2 title "{/Symbol L} = 450 MeV", \
	"./data/ppn-CUTNUM500-40.0-MeV_N4LO-incl-0deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 3 title "{/Symbol L} = 500 MeV",\
	"./data/ppn-CUTNUM550-40.0-MeV_N4LO-incl-0deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 4 title "{/Symbol L} = 500 MeV"

# unset ylabel
# unset ytics

# set format y ""

set title "{{/Symbol q} = 60 deg}" font ",6"  offset 0,-0.5,0
plot "./data/ppn-CUTNUM400-40.0-MeV_N4LO-incl-60deg-2020-01-22" every ::9 u 2:($3*10000)  w l ls 1 title "{/Symbol L} = 400 MeV",\
	"./data/ppn-CUTNUM450-40.0-MeV_N4LO-incl-60deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 2 title "{/Symbol L} = 450 MeV", \
	"./data/ppn-CUTNUM500-40.0-MeV_N4LO-incl-60deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 3 title "{/Symbol L} = 500 MeV",\
	"./data/ppn-CUTNUM550-40.0-MeV_N4LO-incl-60deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 4 title "{/Symbol L} = 500 MeV"

set title "{{/Symbol q} = 120 deg}" font ",6"  offset 0,-0.5,0
plot "./data/ppn-CUTNUM400-40.0-MeV_N4LO-incl-120deg-2020-01-22" every ::9 u 2:($3*10000)  w l ls 1 title "{/Symbol L} = 400 MeV",\
	"./data/ppn-CUTNUM450-40.0-MeV_N4LO-incl-120deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 2 title "{/Symbol L} = 450 MeV", \
	"./data/ppn-CUTNUM500-40.0-MeV_N4LO-incl-120deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 3 title "{/Symbol L} = 500 MeV",\
	"./data/ppn-CUTNUM550-40.0-MeV_N4LO-incl-120deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 4 title "{/Symbol L} = 500 MeV"

set title "{{/Symbol q} = 180 deg}" font ",6"  offset 0,-0.5,0
plot "./data/ppn-CUTNUM400-40.0-MeV_N4LO-incl-180deg-2020-01-22" every ::9 u 2:($3*10000)  w l ls 1 title "{/Symbol L} = 400 MeV",\
	"./data/ppn-CUTNUM450-40.0-MeV_N4LO-incl-180deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 2 title "{/Symbol L} = 450 MeV", \
	"./data/ppn-CUTNUM500-40.0-MeV_N4LO-incl-180deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 3 title "{/Symbol L} = 500 MeV",\
	"./data/ppn-CUTNUM550-40.0-MeV_N4LO-incl-180deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 4 title "{/Symbol L} = 500 MeV"




#################E = 120 MeV

set key font ",6"
set title "{{/Symbol q} = 0 deg}" font ",6"  offset 0,-0.5,0
plot "./data/ppn-CUTNUM400-120.0-MeV_N4LO-incl-0deg-2020-01-22" every ::9 u 2:($3*10000)  w l ls 1 title "{/Symbol L} = 40 MeV",\
	"./data/ppn-CUTNUM450-120.0-MeV_N4LO-incl-0deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 2 title "{/Symbol L} = 450 MeV", \
	"./data/ppn-CUTNUM500-120.0-MeV_N4LO-incl-0deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 3 title "{/Symbol L} = 500 MeV",\
	"./data/ppn-CUTNUM550-120.0-MeV_N4LO-incl-0deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 4 title "{/Symbol L} = 500 MeV"

# unset ylabel
# unset ytics

# set format y ""
set nokey

set title "{{/Symbol q} = 60 deg}" font ",6"  offset 0,-0.5,0
plot "./data/ppn-CUTNUM400-120.0-MeV_N4LO-incl-60deg-2020-01-22" every ::9 u 2:($3*10000)  w l ls 1 title "{/Symbol L} = 400 MeV",\
	"./data/ppn-CUTNUM450-120.0-MeV_N4LO-incl-60deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 2 title "{/Symbol L} = 450 MeV", \
	"./data/ppn-CUTNUM500-120.0-MeV_N4LO-incl-60deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 3 title "{/Symbol L} = 500 MeV",\
	"./data/ppn-CUTNUM550-120.0-MeV_N4LO-incl-60deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 4 title "{/Symbol L} = 500 MeV"

set title "{{/Symbol q} = 120 deg}" font ",6"  offset 0,-0.5,0
plot "./data/ppn-CUTNUM400-120.0-MeV_N4LO-incl-120deg-2020-01-22" every ::9 u 2:($3*10000)  w l ls 1 title "{/Symbol L} = 400 MeV",\
	"./data/ppn-CUTNUM450-120.0-MeV_N4LO-incl-120deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 2 title "{/Symbol L} = 450 MeV", \
	"./data/ppn-CUTNUM500-120.0-MeV_N4LO-incl-120deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 3 title "{/Symbol L} = 500 MeV",\
	"./data/ppn-CUTNUM550-120.0-MeV_N4LO-incl-120deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 4 title "{/Symbol L} = 500 MeV"

set title "{{/Symbol q} = 180 deg}" font ",6"  offset 0,-0.5,0
plot "./data/ppn-CUTNUM400-120.0-MeV_N4LO-incl-180deg-2020-01-22" every ::9 u 2:($3*10000)  w l ls 1 title "{/Symbol L} = 400 MeV",\
	"./data/ppn-CUTNUM450-120.0-MeV_N4LO-incl-180deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 2 title "{/Symbol L} = 450 MeV", \
	"./data/ppn-CUTNUM500-120.0-MeV_N4LO-incl-180deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 3 title "{/Symbol L} = 500 MeV",\
	"./data/ppn-CUTNUM550-120.0-MeV_N4LO-incl-180deg-2020-01-22" every ::9 u 2:($3*10000)  w line ls 4 title "{/Symbol L} = 500 MeV"

unset multiplot