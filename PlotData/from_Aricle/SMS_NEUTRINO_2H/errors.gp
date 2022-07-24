set term pdf enhanced size 5.00in, 2.70in
  set style line 1 lt rgb "orange" lw 2 dt 1 pt 13 pi -10 ps 0.3 ## LO
  set style line 2 lt rgb "blue" lw 2 dt 1  pt 7 pi -10 ps 0.3 ## NLO
  set style line 3 lc rgb "green" lw 2 dt 2 ## N2LO
  set style line 4 lc rgb "red" lw 2 dt 3 ## N3LO
  set style line 5 lc rgb "black" lw 2 dt 1 ## N4LO
  set style line 6 lc rgb "cyan" lw 2 dt 5  ## N4LO+
  set style line 7 lc rgb "violet" lw 2  dt 4 ## AV18
  set pointintervalbox 0.02  ## interval to a point

set output "errors.pdf"
set multiplot
	set lmargin at screen 0.12
	set rmargin at screen 0.535
	set bmargin at screen 0.20
	set tmargin at screen 0.95

set size  0.5,1
set origin 0,0


set key top right
set key Left
set xrange [0:200]
set xtics 0,30, 199
set mxtics 8
set xlabel "E [MeV]"
set ylabel "{/Symbol d}{/Symbol s}_{tot} [%]" offset 1.5,0
# set format y '1e+{%T}'
plot \
	"./data/errors.dat" u 1:2 title "{/Symbol n}_e -> ~{/Symbol n}{.5-}_e" w l ls 2, \
	"./data/errors.dat" u 1:3 title "~{/Symbol n}{.5-}_e -> ~{/Symbol n}{.5-}_e" w l ls 4 dt 4, \
	"./data/errors.dat" u 1:4 title "~{/Symbol n}{.5-}_e -> e^+" w l ls 3

	set lmargin at screen 0.535
	set rmargin at screen 0.92
set xtics 0,30,200
# set format y ""
unset ytics	
set y2tics
unset ylabel
set key bottom left


plot "./data/errors.dat" u 1:5 w l ls 2 title "{/Symbol n}_e -> ~{/Symbol n}{.5-}_e",\
		"./data/errors.dat" u 1:6 w l ls 4 dt 4 title "~{/Symbol n}{.5-}_e -> ~{/Symbol n}{.5-}_e",\
		"./data/errors.dat" u 1:7 w l ls 3 title "~{/Symbol n}{.5-}_e -> e^+"
