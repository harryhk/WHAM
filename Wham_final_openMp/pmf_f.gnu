set term postscript enhanced solid color eps 
set output "pmf_f40000.eps"

set xlabel "z nm"
set ylabel "PMF KJ / mol "

plot 'pmf_f40000_3.xvg' u 1:2 w lp, 'pmf_f40000_2.xvg' u 1:2 w lp, 'pmf_f40000_1.xvg' u 1:2 w lp, 0


