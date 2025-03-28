stats 'select2/etb.csv' index 3 using 2 name 'VBE' nooutput

set grid lw 1.5
set xtics ("W" 0, "L" 20, "G" 45, "X" 74, "U" 84, "K" 85, "G" 116)
set ylabel "Energy (eV)"
set xlabel "k-path"
set key top right

plot 'dft.csv' using 2 with lines lw 1.5 title 'DFT', \
     'select2/etb.csv' using ($2-VBE_max) with lines lw 1.5 title 'ETB'

reset
