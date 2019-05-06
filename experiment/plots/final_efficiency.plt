set terminal postscript eps enhanced color "Helvetica, 22"
set size ratio 0.5
set key top right font "Helvetica, 18"
set title font "Helvetica, 26"

#Griglia:
set style line 1 lc rgb '#808080' lt 0 lw 2
# set grid back ls 1
set ylabel "Efficiency"

load 'color.pal'

set style increment default
set style histogram clustered gap 1 title textcolor lt -1
# set datafile missing '-'
set style data histograms
# set xtics border in scale 0,0 nomirror rotate by -45  autojustify
# set xtics  norangelimit 
# set xtics   ()
set title "Best efficient matcher results for each matcher" 
set xrange [ * : * ] noreverse writeback
set x2range [ * : * ] noreverse writeback
set yrange [ 0.00 : 50. ] noreverse writeback
set y2range [ * : * ] noreverse writeback
set zrange [ * : * ] noreverse writeback
set cbrange [ * : * ] noreverse writeback
set rrange [ * : * ] noreverse writeback

set output "../final_efficiency.eps"

plot '../data/final_results.dat' using ($2/($5/1000/60)):xtic(1) t "sequential", '' u ($2/($6/1000/60)) t "parallel"
