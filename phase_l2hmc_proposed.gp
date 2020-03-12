set palette model HSV
set palette rgb 3,2,2
set cbrange [-pi:pi]
set xrange [-pi:pi]
flip(x) = x>pi ? x-2.*pi : x
set title 'L2HMC proposed points'
set xlabel '$x$'
set ylabel '$p$'
set cblabel 'initial $x$'
unset key
plot 'ergo_l2hmc.out' u (flip($6)):7:(0.4*cos($3*0.25)):(flip($2)) w p pt 7 ps var lc pal z
