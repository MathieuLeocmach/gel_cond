#!/bin/bash

[ $# -eq 0 ] && { echo "Usage: $0 size" ; exit 1; }

mkdir $1

#run povray on percolation sequence
for i in $(seq 8 120); do 
    povray +O$1/$(printf "155C_p2size_t%03d.png" $i) +W$1 +H$1 -D $(printf "155C_percolation_1645_p2size_t%03d.pov" $i); 
done

#run povray on ageing sequence, so that numbers follow percolation sequence
for i in $(seq 0 20); do 
    povray +O$1/$(printf "155C_p2size_t%03d.png" $(($i+121))) +W$1 +H$1 -D $(printf "155C_1715_ageing_p2size_t%02d.pov" $i); 
done

#run ffmpeg to make the video
cd $1
ffmpeg -framerate 24 -pattern_type glob -i "155C_*.png" -vframes 130 -vcodec libx264 -preset veryslow -crf 24 -threads 0 output$1.mp4
cd ..
