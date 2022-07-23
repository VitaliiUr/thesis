set term pdf enhanced size 5.00in, 2.70in
	
  set style line 1 lt rgb "orange" lw 2 dt 1 pt 13 pi -10 ps 0.3 ## LO
  set style line 2 lt rgb "blue" lw 2 dt 1  pt 7 pi -10 ps 0.3 ## NLO
  set style line 3 lc rgb "green" lw 2 dt 2 ## N2LO
  set style line 4 lc rgb "red" lw 2 dt 3 ## N3LO
  set style line 5 lc rgb "black" lw 2 dt 1 ## N4LO
  set style line 6 lc rgb "cyan" lw 2 dt 5  ## N4LO+
  set style line 7 lc rgb "violet" lw 2  dt 4 ## AV18
  set pointintervalbox 0.02  ## interval to a point

set output "Sigmatot_CC_antineutrino_all.pdf"
set multiplot 
	set lmargin at screen 0.14
	set rmargin at screen 0.55
	set bmargin at screen 0.20
	set tmargin at screen 0.95

set size  0.5,1
set origin 0,0

set tics scale 2
# set key at 190,3e8 Left font ",8"
unset key
set xrange [0:200]
set xtics 0,30, 199
set mxtics 3
set xlabel "E_{~{/Symbol n}{.5-}} [MeV]"
set ylabel " {/Symbol s}_{tot} [fm^2] " #offset 2.5,0
set format y '10^{%T}'
set label 1 "(a)" at 20, graph 0.87
set ytics offset 0.5,0

set logscale y
set ytics 1e4, 10, 1e9


plot \
	"./data/antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.0.2.txt" u 1:3 title "LO" w lp ls 1, \
	"./data/antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.1.2.txt" u 1:3 title "NLO" w lp ls 2, \
	"./data/antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.2.2.txt" u 1:3 title "N2LO" w l ls 3, \
	"./data/antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.3.2.txt" u 1:3 title "N3LO" w l ls 4, \
	"./data/antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.4.2.txt" u 1:3 title "N4LO" w l ls 5,\
	"./data/antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.5.2.txt" u 1:3 title "N4LO+" w l ls 6,\
	"./data/antineutrino_CC_deuteron_xs_results_2019_from_table.txt"         u 1:3 title "AV18" w l ls 7
	
	set lmargin at screen 0.55
	set rmargin at screen 0.95
set xtics 0,30,200
set format y ""	
unset ylabel

set label 1 "(b)"
plot "./data/antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.4.1.txt" u 1:3 w l ls 6 title "400 MeV",\
	"./data/antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.4.2.txt" u 1:3 w l ls 5 title "450 MeV",\
	"./data/antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.4.3.txt" u 1:3 w l ls 3 title "500 MeV",\
	"./data/antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.4.4.txt" u 1:3 w l ls 4 title "50 MeV",\
	"./data/antineutrino_CC_deuteron_xs_results_2019_from_table.txt" u 1:3 w l ls 7 title "AV18"
unset label 1
set nokey
unset ylabel
unset xlabel
set lmargin 0.3
set rmargin 0.7
set bmargin 0.2
set tmargin 0.5
set size 0.27, 0.37
set origin 0.267, 0.275
# clear

set autoscale
unset logscale y
# set xrange [120:130]
set xrange [90:95.99]
set yrange [2.2e8:2.41e8]
set xtics 89,2, 97 offset 0,0.5
# set xtics 2 font ",8" offset 0,0.5
set format y "%.1t{x}10^{%T}"
set ytics 0.1e8 font ",11" offset 0.6,0
set mytics 2
set mxtics 4

plot \
	"./data/antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.0.2.txt" u 1:3 title "LO" w lp ls 1 pi -1, \
	"./data/antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.1.2.txt" u 1:3 title "NLO" w lp ls 2 pi -1, \
	"./data/antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.2.2.txt" u 1:3 title "N2LO" w l ls 3, \
	"./data/antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.3.2.txt" u 1:3 title "N3LO" w l ls 4, \
	"./data/antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.4.2.txt" u 1:3 title "N4LO" w l ls 5,\
	"./data/antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.5.2.txt" u 1:3 title "N4LO+" w l ls 6,\
	"./data/antineutrino_CC_deuteron_xs_results_2019_from_table.txt"         u 1:3 title "AV18" w l ls 7



set lmargin 0.3
set rmargin 0.7
set bmargin 0.2
set tmargin 0.5
set size 0.27, 0.37
set origin 0.667, 0.275

plot "./data/antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.4.1.txt" u 1:3 w l ls 6 title "400 MeV",\
	"./data/antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.4.2.txt" u 1:3 w l ls 5 title "450 MeV",\
	"./data/antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.4.3.txt" u 1:3 w l ls 3 title "500 MeV",\
	"./data/antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.4.4.txt" u 1:3 w l ls 4 title "50 MeV",\
	"./data/antineutrino_CC_deuteron_xs_results_2019_from_table.txt" u 1:3 w l ls 7 title "AV18"
