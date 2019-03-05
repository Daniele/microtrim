#!/usr/bin/gnuplot -persist

set terminal pdfcairo enhanced font 'Verdana,10'
set output 'v3_parallel-scalability.pdf'

# set title "Latency Benchmark (100 pipelines, 8 actor each, 100 big-actors, vector 5000).\n"

set xlabel "N°. of workers"
set ylabel "Scalability"
set xrange[1:23]
set xtics 1

stats 'v3_parallel.dat' every ::1::23 u 5 nooutput name 'Tc_'
set yrange[1:23]

set grid

ideal(x) = x
plot 'v3_parallel.dat' every ::1::23 u 2:(Tc_max/$5) w linespoints t 'V3 chunk1000',\
      ideal(x) title 'ideal' with lines

# pause -1
