set palette model HSV
set palette rgb 3,2,2
set cbrange [-pi:pi]
set xrange [-pi:pi]
flip(x) = x>pi ? x-2.*pi : x
set title 'HMC proposed points'
set xlabel '$x$'
set ylabel '$p$'
set cblabel 'initial $x$'
unset key
plot 'ergo_hmc.out' u (flip($5)):6:(0.4*cos($3*0.25)):(flip($2)) w p pt 7 ps var lc pal z
