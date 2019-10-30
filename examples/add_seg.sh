#!/usr/bin/env expect

for {set x 1} {$x<5} {incr x} {
    eval spawn "panopticon book cnv-wizard"
    expect "segmentation data"
    send "./example_data/example.loom\n"
    expect "add to loom file"
    send "./example_data/Example_Seg_$x.seg\n"
    expect "Which column"
    send "Chromosome\n"
    expect "Which column"
    send "Start.bp\n"
    expect "Which column"
    send "End.bp\n"
    expect "Which column"
    send "tau\n"
    expect "Was that"
    send "y\n"
    expect "like to label"
    send "Seg$x\n"
    expect "complete"

}
