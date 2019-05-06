set terminal postscript eps enhanced color "Helvetica, 22"
set key top right font "Helvetica, 18"
set title font "Helvetica, 26"

set output '../accuracy_detail.eps'
set boxwidth 0.9 absolute
set style fill   solid 1.00 border lt -1
set key fixed right top vertical Right noreverse noenhanced autotitle nobox
set style increment default
set style histogram clustered gap 1 title textcolor lt -1
set datafile missing '-'
set style data histograms
set xtics border in scale 0,0 nomirror rotate by -45  autojustify
set xtics  norangelimit 
set xtics   ()
set title "HTSeq results for each matcher" 
set xrange [ * : * ] noreverse writeback
set x2range [ * : * ] noreverse writeback
set yrange [ 32. : 34.5 ] noreverse writeback
set y2range [ * : * ] noreverse writeback
set zrange [ * : * ] noreverse writeback
set cbrange [ * : * ] noreverse writeback
set rrange [ * : * ] noreverse writeback
set ylabel "Percentage %"
## Last datafile plotted: "immigration.dat"
plot '../data/accuracy_detail.dat' using 2:xtic(1) ti col
