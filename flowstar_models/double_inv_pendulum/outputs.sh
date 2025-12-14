#!/bin/bash

cd outputs
gnuplot cart_motion.plt
gnuplot lqr_angles.plt

convert cart_motion.eps cart_motion.png
convert lqr_angles.eps lqr_angles.png

cd ..
