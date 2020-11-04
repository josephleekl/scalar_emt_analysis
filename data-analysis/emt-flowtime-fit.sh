#!/usr/bin/env bash 
#
###############################################
#
# Perform fit of c3 against flowtime for each 
# ensemble and save plots
#
###############################################



emt-flowtime-c3-fit data/param/flowtime-fit/flowtime-c3-fit-par-g0.1-L64-m2-0.0305.xml -d --save-plot 
emt-flowtime-c3-fit data/param/flowtime-fit/flowtime-c3-fit-par-g0.1-L64-m2-0.031.xml -d --save-plot
emt-flowtime-c3-fit data/param/flowtime-fit/flowtime-c3-fit-par-g0.1-L128-m2-0.0305.xml -d --save-plot
emt-flowtime-c3-fit data/param/flowtime-fit/flowtime-c3-fit-par-g0.1-L128-m2-0.031.xml -d --save-plot
emt-flowtime-c3-fit data/param/flowtime-fit/flowtime-c3-fit-par-g0.1-L256-m2-0.0305.xml -d --save-plot
emt-flowtime-c3-fit data/param/flowtime-fit/flowtime-c3-fit-par-g0.1-L256-m2-0.031.xml -d --save-plot

emt-flowtime-c3-fit data/param/flowtime-fit/flowtime-c3-fit-par-g0.2-L64-m2-0.061.xml -d --save-plot
emt-flowtime-c3-fit data/param/flowtime-fit/flowtime-c3-fit-par-g0.2-L64-m2-0.062.xml -d --save-plot
emt-flowtime-c3-fit data/param/flowtime-fit/flowtime-c3-fit-par-g0.2-L128-m2-0.061.xml -d --save-plot
emt-flowtime-c3-fit data/param/flowtime-fit/flowtime-c3-fit-par-g0.2-L128-m2-0.062.xml -d --save-plot
emt-flowtime-c3-fit data/param/flowtime-fit/flowtime-c3-fit-par-g0.2-L256-m2-0.061.xml -d --save-plot
emt-flowtime-c3-fit data/param/flowtime-fit/flowtime-c3-fit-par-g0.2-L256-m2-0.062.xml -d --save-plot



emt-flowtime-c3-fit data/param/flowtime-fit/flowtime-c3-fit-par-g0.3-L64-m2-0.091.xml -d --save-plot
emt-flowtime-c3-fit data/param/flowtime-fit/flowtime-c3-fit-par-g0.3-L64-m2-0.092.xml -d --save-plot
emt-flowtime-c3-fit data/param/flowtime-fit/flowtime-c3-fit-par-g0.3-L128-m2-0.091.xml -d --save-plot
emt-flowtime-c3-fit data/param/flowtime-fit/flowtime-c3-fit-par-g0.3-L128-m2-0.092.xml -d --save-plot
emt-flowtime-c3-fit data/param/flowtime-fit/flowtime-c3-fit-par-g0.3-L256-m2-0.091.xml -d --save-plot
emt-flowtime-c3-fit data/param/flowtime-fit/flowtime-c3-fit-par-g0.3-L256-m2-0.092.xml -d --save-plot