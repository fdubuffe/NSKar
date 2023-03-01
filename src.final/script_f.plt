#set terminal postscript eps enhanced color
set terminal png
set output "input.png"
set size ratio -1
set origin 0,0
set xrange [0:1]
set yrange [0:1]
#set cbrange [-3.1368634184076867:+3.1368634184076867]
#set cbrange [vmin:vmax]
set palette defined (0. "cyan", 0.25 "blue", 0.5 "black", 0.75 "red", 1 "yellow")
set colorbox
set pm3d map
set pm3d flush begin ftriangles scansforward interpolate 10,10
splot 'input.plt' w pm3d notitle
