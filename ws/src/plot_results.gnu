set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output 'stability_analysis.png'

set multiplot layout 2,2 title 'Stability Analysis (U=4.47214 m/s)'

set title 'Depth and Pitch vs Time'
set xlabel 'Time (s)'
set ylabel 'Depth z (m)'
set y2label 'Pitch θ (rad)'
set y2tics
plot 'trajectory_data.txt' using 1:3 with lines lw 2 title 'Depth (z)', \
     'trajectory_data.txt' using 1:4 with lines lw 2 axes x1y2 title 'Pitch (θ)'

set title 'Phase Plane: Depth'
set xlabel 'Depth z (m)'
set ylabel 'Vertical velocity ż (m/s)'
unset y2tics
plot 'trajectory_data.txt' using 3:5 with lines lw 2 title 'Phase trajectory'

set title 'Phase Plane: Pitch Angle'
set xlabel 'Pitch θ (rad)'
set ylabel 'Angular velocity θ̇ (rad/s)'
plot 'trajectory_data.txt' using 4:6 with lines lw 2 title 'Phase trajectory'

set title 'Spatial Trajectory'
set xlabel 'X position (m)'
set ylabel 'Depth z (m)'
plot 'trajectory_data.txt' using 2:3 with lines lw 2 title 'Trajectory'

unset multiplot
set output
