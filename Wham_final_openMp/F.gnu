set term postscript enhanced solid color eps 
set output "pmf_z4000.eps"

set xlabel "z nm"
set ylabel "PMF KJ / mol "

plot 'pmf_z40000_3.xvg' u 1:2 w lp, 'pmf_z40000_2.xvg' u 1:2 w lp, 'pmf_z40000_1.xvg' u 1:2 w lp, 0


