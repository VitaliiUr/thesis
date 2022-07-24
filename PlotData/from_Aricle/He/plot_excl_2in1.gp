set term pdf enhanced size 5.00in, 2.70in


  set style line 1 lt rgb "orange" lw 2 dt 1 pt 13 pi -10 ps 0.3 ## LO
  set style line 2 lt rgb "blue" lw 2 dt 1  pt 7 pi -10 ps 0.3 ## NLO
  set style line 3 lc rgb "green" lw 2 dt 2 ## N2LO
  set style line 4 lc rgb "red" lw 2 dt 3 ## N3LO
  set style line 5 lc rgb "black" lw 2 dt 1 ## N4LO
  set style line 6 lc rgb "cyan" lw 2 dt 5  ## N4LO+
  set style line 7 lc rgb "violet" lw 2  dt 4 ## AV18
  set pointintervalbox 0.02  ## interval to a point


MP_LEFT = .14
MP_RIGHT = .95
MP_BOTTOM = .16
MP_TOP = .95
MP_GAP = 0.00005


set output "Excl120MeV_15_0_15_180_cutnum450_all.pdf"
set multiplot layout 1,2 margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP

set tics scale 2
set format y "%g"
# set key center top
set key at graph 0.75, graph 0.9 # font ",14"
# unset key
set ylabel "d^5{/Symbol s}/d{/Symbol W}_1d{/Symbol W}_2dS [{/Symbol m}b sr^{-2} MeV^{-1}]" #offset 2.5,0
# set title "E = 120 MeV, {/Symbol Q}1_{lab} = 15^0, {/Symbol f}1_{lab} = 0^0, {/Symbol Q}2_{lab} = 15^0, {/Symbol f}2_{lab} = 180^0"
set xlabel "S [MeV]"
set xrange [0:120]
set xtics 0,20,119
set yrange [0:0.025]
set mxtics 2
set mytics 2
set label 1 "(a)" at 13, graph 0.12

plot "./data/ppn-CUTNUM450-120.0-MeV_LO" u 2:($3*10000) w lp ls 1 title "LO",\
		"./data/ppn-CUTNUM450-120.0-MeV_NLO" u 2:($3*10000) w lp ls 2 title "NLO",\
		"./data/ppn-CUTNUM450-120.0-MeV_N2LO" u 2:($3*10000) w l ls 3 title "N^2LO",\
		"./data/ppn-CUTNUM450-120.0-MeV_N3LO" u 2:($3*10000) w l ls 4 title "N^3LO",\
		"./data/ppn-CUTNUM450-120.0-MeV_N4LO" u 2:($3*10000) w l ls 5 title "N^4LO",\
		"./data/ppn-CUTNUM450-120.0-MeV_N4LO+" u 2:($3*10000) w l ls 6 title "N^4LO+"

set xtics 0,20,120
set format y ""	
# set key center top
unset ylabel
# set title "E = 120 MeV, {/Symbol Q}1_{lab} = 15^0, {/Symbol f}1_{lab} = 0^0, {/Symbol Q}2_{lab} = 15^0, {/Symbol f}2_{lab} = 180^0"
# set xlabel "E_p"
set label 1 "(b)"

plot "./data/ppn-CUTNUM400-120.0-MeV_N4LO" u 2:($3*10000)  w l ls 6 title "400 MeV",\
	"./data/ppn-CUTNUM450-120.0-MeV_N4LO" u 2:($3*10000)  w l ls 5 title "450 MeV", \
	"./data/ppn-CUTNUM500-120.0-MeV_N4LO" u 2:($3*10000)  w l ls 3 title "500 MeV",\
	"./data/ppn-CUTNUM550-120.0-MeV_N4LO" u 2:($3*10000)  w l ls 4 title "550 MeV"

unset multiplot

set output "Excl40MeV_15_0_15_180_cutnum450_all.pdf"
set multiplot layout 1,2 margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP

set format y "%g"
# set key at 16,0.0295
# unset key
set ylabel "d^5{/Symbol s}/d{/Symbol W}_1d{/Symbol W}_2dS [{/Symbol m}b sr^{-2} MeV^{-1}]" #offset 2.5,0
# set title "E = 120 MeV, {/Symbol Q}1_{lab} = 15^0, {/Symbol f}1_{lab} = 0^0, {/Symbol Q}2_{lab} = 15^0, {/Symbol f}2_{lab} = 180^0"
set xlabel "S [MeV]"
set autoscale
set xtics autofreq
set xrange [0:30]
set yrange [0:0.026]
set xtics 0,5,29.9
set label 1 "(a)" at 5, graph 0.12

plot "./data/ppn-CUTNUM450-40.0-MeV_LO" u 2:($3*10000) w lp ls 1 title "LO",\
		"./data/ppn-CUTNUM450-40.0-MeV_NLO" u 2:($3*10000) w lp ls 2 title "NLO",\
		"./data/ppn-CUTNUM450-40.0-MeV_N2LO" u 2:($3*10000) w l ls 3 title "N^2LO",\
		"./data/ppn-CUTNUM450-40.0-MeV_N3LO" u 2:($3*10000) w l ls 4 title "N^3LO",\
		"./data/ppn-CUTNUM450-40.0-MeV_N4LO" u 2:($3*10000) w l ls 5 title "N^4LO",\
		"./data/ppn-CUTNUM450-40.0-MeV_N4LO+" u 2:($3*10000) w l ls 6 title "N^4LO+"

set xtics 0,5,30
set format y ""	
# set key at 16,0.029
unset ylabel
# set title "E = 120 MeV, {/Symbol Q}1_{lab} = 15^0, {/Symbol f}1_{lab} = 0^0, {/Symbol Q}2_{lab} = 15^0, {/Symbol f}2_{lab} = 180^0"
# set xlabel "E_p"
set label 1 "(b)"

plot "./data/ppn-CUTNUM400-40.0-MeV_N4LO" u 2:($3*10000)  w l ls 6 title "400 MeV",\
	"./data/ppn-CUTNUM450-40.0-MeV_N4LO" u 2:($3*10000)  w l ls 5 title "450 MeV", \
	"./data/ppn-CUTNUM500-40.0-MeV_N4LO" u 2:($3*10000)  w l ls 3 title "500 MeV",\
	"./data/ppn-CUTNUM550-40.0-MeV_N4LO" u 2:($3*10000)  w l ls 4 title "550 MeV"