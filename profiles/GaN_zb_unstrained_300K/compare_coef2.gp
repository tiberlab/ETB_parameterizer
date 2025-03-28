# Set up the layout with extra space for the legend
#set lmargin 10
#set rmargin 12
# Set up the multiplot layout
set multiplot layout 4,2

# Common settings
set grid
set yrange [-0.01:1.01]
set xtics ("W" 0, "L" 20, "G" 45, "X" 74, "U" 84, "K" 85, "G" 116)
#set ylabel "Orbital character"

#set size 1.2,1  # Make the plot 20% wider to accommodate legend
unset key
#set key at screen 1,0.8 vertical box spacing 1.2  # Position legend in the extra space


set title "Lowest CB, DFT"
plot 'dft.csv' index 4 using 3 with lines lc rgb "red" title "s-Ga", \
    '' index 4 using 4 with lines lc rgb "black" title "p-Ga", \
    '' index 4 using 5 with lines lc rgb "white" title "d-Ga", \
    '' index 4 using 6 with lines dashtype 2 lc rgb "red" title "s-N", \
    '' index 4 using 7 with lines dashtype 2 lc rgb "black" title "p-N", \
    '' index 4 using 8 with lines dashtype 2 lc rgb "white" title "d-N"


set title "Lowest CB, ETB"
plot 'select2/etb.csv' index 4 using 3 with lines lc rgb "red" title "s-Ga", \
    '' index 4 using 4 with lines lc rgb "black" title "p-Ga", \
    '' index 4 using 5 with lines lc rgb "white" title "d-Ga", \
    '' index 4 using 6 with lines dashtype 2 lc rgb "red" title "s-N", \
    '' index 4 using 7 with lines dashtype 2 lc rgb "black" title "p-N", \
    '' index 4 using 8 with lines dashtype 2 lc rgb "white" title "d-N"

set key at screen 0.3, 0.66

set title "HH band, DFT"
plot 'dft.csv' index 3 using 3 with lines lc rgb "red" title "s-Ga", \
    '' index 3 using 4 with lines lc rgb "black" title "p-Ga", \
    '' index 3 using 5 with lines lc rgb "white" title "d-Ga", \
    '' index 3 using 6 with lines dashtype 2 lc rgb "red" title "s-N", \
    '' index 3 using 7 with lines dashtype 2 lc rgb "black" title "p-N", \
    '' index 3 using 8 with lines dashtype 2 lc rgb "white" title "d-N"

unset key 

set title "HH band, ETB"
plot 'select2/etb.csv' index 3 using 3 with lines lc rgb "red" title "s-Ga", \
    '' index 3 using 4 with lines lc rgb "black" title "p-Ga", \
    '' index 3 using 5 with lines lc rgb "white" title "d-Ga", \
    '' index 3 using 6 with lines dashtype 2 lc rgb "red" title "s-N", \
    '' index 3 using 7 with lines dashtype 2 lc rgb "black" title "p-N", \
    '' index 3 using 8 with lines dashtype 2 lc rgb "white" title "d-N"

set title "LH band, DFT"
plot 'dft.csv' index 2 using 3 with lines lc rgb "red" title "s-Ga", \
    '' index 2 using 4 with lines lc rgb "black" title "p-Ga", \
    '' index 2 using 5 with lines lc rgb "white" title "d-Ga", \
    '' index 2 using 6 with lines dashtype 2 lc rgb "red" title "s-N", \
    '' index 2 using 7 with lines dashtype 2 lc rgb "black" title "p-N", \
    '' index 2 using 8 with lines dashtype 2 lc rgb "white" title "d-N"

set title "LH band, ETB"
plot 'select2/etb.csv' index 2 using 3 with lines lc rgb "red" title "s-Ga", \
    '' index 2 using 4 with lines lc rgb "black" title "p-Ga", \
    '' index 2 using 5 with lines lc rgb "white" title "d-Ga", \
    '' index 2 using 6 with lines dashtype 2 lc rgb "red" title "s-N", \
    '' index 2 using 7 with lines dashtype 2 lc rgb "black" title "p-N", \
    '' index 2 using 8 with lines dashtype 2 lc rgb "white" title "d-N"

set title "SO band, DFT"
plot 'dft.csv' index 1 using 3 with lines lc rgb "red" title "s-Ga", \
    '' index 1 using 4 with lines lc rgb "black" title "p-Ga", \
    '' index 1 using 5 with lines lc rgb "white" title "d-Ga", \
    '' index 1 using 6 with lines dashtype 2 lc rgb "red" title "s-N", \
    '' index 1 using 7 with lines dashtype 2 lc rgb "black" title "p-N", \
    '' index 1 using 8 with lines dashtype 2 lc rgb "white" title "d-N"

set title "SO band, ETB"
plot 'select2/etb.csv' index 1 using 3 with lines lc rgb "red" title "s-Ga", \
    '' index 1 using 4 with lines lc rgb "black" title "p-Ga", \
    '' index 1 using 5 with lines lc rgb "white" title "d-Ga", \
    '' index 1 using 6 with lines dashtype 2 lc rgb "red" title "s-N", \
    '' index 1 using 7 with lines dashtype 2 lc rgb "black" title "p-N", \
    '' index 1 using 8 with lines dashtype 2 lc rgb "white" title "d-N"

unset multiplot

reset
