set term pdfcairo enhanced size 5.20in, 1.90in

set style line 1 lw 2 ps 0.5 lc rgb "dark-violet" dt 2  #SCS
set style line 2 lw 2 pt 1 ps 0.5 lc rgb "dark-green" dt 1 #SMS


MP_LEFT = .08
MP_RIGHT = .95
MP_BOTTOM = .18
MP_TOP = .9
MP_GAP = 0.06

set output "Diff.pdf"

set multiplot layout 1,3 margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP
		# title "E=30 MeV"
set size square
# set format '$%g$'
# set label 2 "Chiral order" at screen 0.4, screen 0.05 font ",10"	
set ylabel "{/Symbol D}{/Symbol s}[{/Symbol m}b sr^{-1}]"	offset 2, 0 font ",12"

# unset key

# set tics scale 2
set xrange [0.5:5.5]
# set yrange [0:0.35]
# set mytics 2
# set mxtics 2
set ytics font ",12"
set xtics ("LO" 1, "NLO" 2, "N^2LO" 3, "N^3LO" 4,"N^4LO" 5) font ",9"
# set grid

set title "E = 30MeV, {/Symbol q}_p = 60^o" font ",10"
set label 1 "(a)" at 1, graph 0.12

plot "diffDegCMS.dat" every ::2 u ($0+1):1 w lp ls 1 title "SCS",\
	"diffDegSMS.dat" every ::2 u ($0+1):1 w lp ls 2 title "SMS"

unset key
# set yrange [0:0.055]
# set ytics 0,0.01,0.05
# unset ylabel
# unset ytics
# set format y ""
set label 1 "(b)"
unset ylabel
set title "E = 100MeV, {/Symbol q}_p = 15^o" font ",10"
plot "diffDegCMS.dat" every ::2 u ($0+1):2 w lp ls 1 title "SCS",\
	"diffDegSMS.dat" every ::2 u ($0+1):2 w lp ls 2 title "SMS"

set label 1 "(c)"

set title "E = 100MeV, {/Symbol q}_p = 150^o" font ",10"
plot "diffDegCMS.dat" every ::2 u ($0+1):3 w lp ls 1 title "SCS",\
	"diffDegSMS.dat" every ::2 u ($0+1):3 w lp ls 2 title "SMS"


unset multiplot
