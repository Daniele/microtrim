set terminal postscript eps enhanced color "Helvetica, 22"
set size ratio 0.5
set key top right font "Helvetica, 18"
set title font "Helvetica, 26"

#Griglia:
#set style line 1 lc rgb '#808080' lt 0 lw 2
# set grid back ls 1
set style fill   solid 1.00 border lt -1
set ylabel "Percentage %"

load 'color.pal'

set style increment default
set style histogram clustered gap 1 title textcolor lt -1
# set datafile missing '-'
set style data histograms
set xtics border in scale 0,0 nomirror rotate by -45  autojustify
# set xtics  norangelimit 
# set xtics   ()
set key outside
set key center bottom
set title "Best efficient matcher results for each matcher" 
set xrange [ * : * ] noreverse writeback
set x2range [ * : * ] noreverse writeback
set yrange [ 0.00 : 100. ] noreverse writeback
set y2range [ * : * ] noreverse writeback
set zrange [ * : * ] noreverse writeback
set cbrange [ * : * ] noreverse writeback
set rrange [ * : * ] noreverse writeback

set output "../final_efficiency.eps"

#((($2/33.91)*((1-($5/1148005))))*100):xtic(1) t "sequential",

plot '../data/final_results.dat' using (($2/33.91)*100):xtic(1) t "accuracy", '' u ((1-($5/1148005))*100) t "sequential time", '' u ((1-($6/1148005))*100) t "parallel time"
