set term pdfcairo# enhanced

set style line 1 lc rgb "dark-violet" dt 2
set style line 2 lc rgb "dark-green" dt 1

MP_LEFT = .1
MP_RIGHT = .95
MP_BOTTOM = .17
MP_TOP = .95
MP_GAP = 0.1

set output "Delta30MeV.pdf"

set multiplot layout 1,2 margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP
		# title "E=30 MeV"

# set format '$%g$'
set xlabel "{/Symbol q}_p [deg]"	
set ylabel "{/Symbol D}{/Symbol s} / {/Symbol s}_{avrg}"	offset 2, 0

set xrange [0:180]
# set yrange [0:0.35]
set xtics 0,30,180
set mytics 2
set mxtics 2


# set title "Chiral order"

plot "deltaData_old_30MeV.dat" u 1:4 w l ls 1 title "SCS", "deltaData_30MeV.dat" u 1:4 w l ls 2 title "SMS"

unset key
set yrange [0:0.055]
set ytics 0,0.01,0.05
# unset ylabel
# unset ytics
# set format y ""
unset ylabel
# set title "Cut-off"
plot "deltaData_old_30MeV.dat" u 1:5 w l ls 1 title "SCS", "deltaData_30MeV.dat" u 1:5 w l ls 2 title "SMS"


unset multiplot

set output "Delta100MeV.pdf"


set multiplot layout 1,2 margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP
		# title "E=100 MeV"

set autoscale y
set ytics autofreq
set key
set ylabel "{/Symbol D}{/Symbol s} / {/Symbol s}_{avrg}"	offset 2, 0

# set title "Chiral order"
plot "deltaData_old_100MeV.dat" u 1:2 w l ls 1 title "SCS", "deltaData_100MeV.dat" u 1:2 w l ls 2 title "SMS"

unset key
# unset ylabel
# unset ytics
# set format y ""
unset ylabel
# set title "Cut-off"
plot "deltaData_old_100MeV.dat" u 1:5 w l ls 1 title "SCS", "deltaData_100MeV.dat" u 1:5 w l ls 2 title "SMS"

