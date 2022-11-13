#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 16:03:10 2022

@author: dyoung
"""
#!/usr/bin/env python3
#coding=utf8
from pygnuplot import gnuplot

# Illustration of object-oriented interface, you can see we only wrap the
# gnuplot script by g.cmd('...') and it's simple and straitfoward if you
# are familar with Gnuplot.
g = gnuplot.Gnuplot()
g.cmd('set terminal pngcairo font "arial,10" fontscale 1.0 size 600, 400')
g.cmd('set output "simple.png"')
g.cmd('set key fixed left top vertical Right noreverse enhanced autotitle box lt black linewidth 1.000 dashtype solid')
g.cmd('set style increment default')
g.cmd('set samples 50, 50')
g.cmd('set title "Simple Plots" ')
g.cmd('set title  font ",20" norotate')
g.cmd('plot [-10:10] sin(x),atan(x),cos(atan(x))')