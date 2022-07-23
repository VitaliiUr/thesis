set term pdfcairo enhanced size 5.00in, 2.00in

set style line 1 lc rgb "dark-violet" lw 2 dt 2 #SCS
set style line 2 lc rgb "dark-green" lw 2 dt 1 #SMS

MP_LEFT = .07
MP_RIGHT = .95
MP_BOTTOM = .24
MP_TOP = .85
MP_GAP = 0.07

set output "Delta_3in1_30MeV.pdf"

set multiplot layout 1,3 margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP
		# title "E=30 MeV"

# set format '$%g$'
set xlabel "{/Symbol q}_p [deg]"	
# set ylabel "{/Symbol D}{/Symbol s} / {/Symbol s}_{avrg}"	offset 2, 0

set tics scale 2
set xrange [0:180]
# set yrange [0:0.35]
set xtics 60
set mytics 2
set mxtics 6
# unset key

# set title "{/Symbol D}{/Symbol s} / {/Symbol s}_{avrg}"
set title "{/Symbol d}{/Symbol s} (chiral order)[%]"

set label 1 "(a)" at 20, graph 0.12 font ",12"

plot "deltaData_old_30MeV.dat" u 1:(100*$4) w l ls 1 title "SCS", "deltaData_30MeV.dat" u 1:(100*$4) w l ls 2 title "SMS"

unset key
unset ylabel
# unset label 1


set yrange [0:5.5]
set ytics 1
set label 1 "(b)" at 80, graph 0.12
# unset ylabel
# unset ytics
# set format y ""
# set title "{/Symbol D}{/Symbol s} / {/Symbol s}_{avrg({/Symbol L})}"
set title "{/Symbol d}{/Symbol s} ({/Symbol L})[%]"
plot "deltaData_old_30MeV.dat" u 1:(100*$5) w l ls 1 title "SCS", "deltaData_30MeV.dat" u 1:(100*$5) w l ls 2 title "SMS"


set autoscale
set xrange [0:180]
set label 1 "(c)"
set ytics autofreq
set title "{/Symbol D}{/Symbol s}_{(N4LO, N3LO)}[{/Symbol m}b sr^{-1}]"
plot "deltaData_old_30MeV.dat" u 1:7 w l ls 1 title "SCS", "deltaData_30MeV.dat" u 1:7 w l ls 2 title "SMS", 

unset multiplot

set output "Delta_3in1_100MeV.pdf"


set multiplot layout 1,3 margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP
		# title "E=100 MeV"

set autoscale y
set ytics 10
set key
set label 1 "(a)" at 20, graph 0.12
# set title "{/Symbol D}{/Symbol s} / {/Symbol s}_{avrg}"
set title "{/Symbol d}{/Symbol s} (chiral order)[%]"
plot "deltaData_old_100MeV.dat" u 1:(100*$4) w l ls 1 title "SCS", "deltaData_100MeV.dat" u 1:(100*$4) w l ls 2 title "SMS"

unset key

# unset ylabel
set ytics 5
set mytics 2
# set format y ""
unset ylabel
set label 1 "(b)"
# set title "{/Symbol D}{/Symbol s} / {/Symbol s}_{avrg({/Symbol L})}"
set title "{/Symbol d}{/Symbol s} ({/Symbol L})[%]"
plot "deltaData_old_100MeV.dat" u 1:(100*$5) w l ls 1 title "SCS", "deltaData_100MeV.dat" u 1:(100*$5) w l ls 2 title "SMS"


set autoscale
set xrange [0:180]
set ytics autofreq
set mytics 2
set title "{/Symbol D}{/Symbol s}_{(N4LO, N3LO)}[{/Symbol m}b sr^{-1}]"
set label 1 "(c)"
plot "deltaData_old_100MeV.dat" u 1:7 w l ls 1 title "SCS", "deltaData_100MeV.dat" u 1:7 w l ls 2 title "SMS", 





# plot "CROSS_30mev.dat" u 1:5 w l title "SMS N3LO",\
 # "CROSS_30mev.dat" u 1:7 w l title "SMS N4LO",\
 # "CROSS_30mev.dat" u 1:10 w l  title "SMS N4LO+",\
 #  "CROSS_old_30mev.dat" u 1:5 w l title "SCS N3LO",\
 #  "CROSS_old_30mev.dat" u 1:5 w l title "SCS N3LO"
