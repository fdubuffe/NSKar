set terminal postscript eps enhanced color
set output "tp.eps"
set size ratio -1
set origin 0,0
set xrange [0:5.6715856]
set yrange [0:1]
set cbrange [-0.5:+0.5]
set palette defined (0. "cyan", 0.25 "blue", 0.5 "black", 0.75 "red", 1 "yellow")
set colorbox
set pm3d map
set pm3d flush begin ftriangles scansforward interpolate 10,10
#splot "fort.99"
splot 'tp01500000.plt' w pm3d notitle

#set terminal postscript landscape
#set output "dessin.ps"
##set terminal jpeg
##set output "dessin.jpeg"
#set size ratio -1
#set origin 0,0
#set xlabel "x"
#set ylabel "z"
#set xrange [0:5.6715856]
#set yrange [0:1]
#set zrange [-0.5:0.5]
#set pm3d map
#set palette defined (0. "cyan", 0.25 "blue", 0.5 "black", 0.75 "red" , 1 "yellow")
#set colorbox
#set lmargin 0
#set title "Plot"
##set pm3d flush begin ftriangles scansforward interpolate 10,10
#splot "fort.99"
