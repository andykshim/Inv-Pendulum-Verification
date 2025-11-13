#!/bin/bash

cd outputs
gnuplot lqr_cart_motion.plt
gnuplot lqr_config_space.plt
gnuplot lqr_phase_portrait.plt

convert lqr_cart_motion.eps lqr_cart_motion.png
convert lqr_config_space.eps lqr_config_space.png
convert lqr_phase_portrait.eps lqr_phase_portrait.png

cd ..
