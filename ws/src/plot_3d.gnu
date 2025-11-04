set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output '3d_trajectory.png'

set title '3D Trajectory: X-Y-Time'
set xlabel 'X position'
set ylabel 'Y position'
set zlabel 'Time'
splot 'trajectory_data.txt' using 2:3:1 with lines lw 2 title 'Trajectory'
