#to be launched on script/ folder
#Plot speedup (i.e. considering the serial execution time)
#Questi sono i grafici da inserire nell'articolo

set terminal postscript eps enhanced color "Helvetica, 22"
set size ratio 0.5
set key top right font "Helvetica, 18"
set title font "Helvetica, 26"

#Griglia:
set style line 1 lc rgb '#808080' lt 0 lw 2
set grid back ls 1
set ylabel "Scalability"
set xlabel "Parallelism Degree"
load 'color.pal'


set title "Scalability"
set output "../final_scal.eps"
set xtics 4
set xrange[1:48]
set ytics 4
set yrange[1:48]
plot\
  "../data/final_ssw.dat"         u 1:(312403/$3)  with lp ls 5 lw 3.5 pt 2 ps 2 title "Striped Smith-Waterman",\
  "../data/final_ndleven.dat"     u 1:(2223405/$3) with lp ls 4 lw 3.5 pt 2 ps 2 title "Damerau-Levenshtein",\
  "../data/final_leven.dat"       u 1:(810127/$3)  with lp ls 3 lw 3.5 pt 2 ps 2 title "Levenshtein",\
  "../data/final_adagen-fast.dat" u 1:(252162/$3)  with lp ls 1 lw 3.5 pt 2 ps 2 title "Adagen fast",\
  "../data/final_adagen.dat"      u 1:(1076678/$3) with lp ls 2 lw 3.5 pt 2 ps 2 title "Adagen",\
