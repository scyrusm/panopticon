#!/usr/bin/env expect

set timeout 120
eval spawn panopticon windowed-mean-expression-clustering --patient 43 --cell_type tumor --complexity_cutoff 1000 --figure_output example.pdf example_data/example.loom
expect "save these clusters"
send "n\n"



