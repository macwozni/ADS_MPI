#!/usr/bin/env python

import sys
from subprocess import Popen, PIPE


text = '''
    set pm3d map
    set pm3d interpolate %(interpolation)s
    set title "%(title)s"
    set terminal png size 900,900
    set output "%(output)s"
    splot "%(file)s"
    set xrange [%(xrange)s]
    set yrange [%(yrange)s]
'''

params = {
    'interpolation': '10,10',
    'xrange': '0:9',
    'yrange': '0:9',
}


name = sys.argv[1]
out = sys.argv[2]
title = 'Random title'

conf = {
    'file': name,
    'output': out,
    'title': title,
}


cmd = text % dict(params, **conf)

p = Popen(['gnuplot'], shell=True, stdin=PIPE)
p.communicate(cmd)

