set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output 'trajectory_plots.png'

set multiplot layout 2,2 title 'Fish Motion Analysis'

set title 'Trajectory X-Y'
set xlabel 'X position'
set ylabel 'Y position'
plot 'trajectory_data.txt' using 2:3 with lines lw 2 title 'Trajectory', \
     'trajectory_data.txt' every 50 using 2:3 with points pt 7 ps 0.5 title 'Points'

set title 'Position vs Time'
set xlabel 'Time'
set ylabel 'Position'
plot 'trajectory_data.txt' using 1:2 with lines lw 2 title 'X position', \
     'trajectory_data.txt' using 1:3 with lines lw 2 title 'Y position'

set title 'Angle vs Time'
set xlabel 'Time'
set ylabel 'Angle (rad)'
plot 'trajectory_data.txt' using 1:4 with lines lw 2 title 'Theta'

set title 'Velocity vs Time'
set xlabel 'Time'
set ylabel 'Velocity'
plot 'trajectory_data.txt' using 1:5 with lines lw 2 title 'Vx', \
     'trajectory_data.txt' using 1:6 with lines lw 2 title 'Vy'

unset multiplot
set output
