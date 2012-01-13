#!/bin/bash

if [ $# -ne 2 ]
then
    echo
    echo "This scripts starts the desired number of stanford parser instances."
    echo "Usage:"
    echo "$0 NUMBER_OF_CORES STARTING_PORT_NUMBER"
    echo
else
    function on_die()
    {  
        echo 'Killing servers...'
        kill $(jobs -pr)
        wait
        exit 0
    }
    
    trap 'on_die' EXIT

    cd ../../3rdparty/stanford-corenlp/
    for p in $(seq $2 1 $(($1+$2-1)))
    do
        python corenlp.py -H 0.0.0.0 -p $p &
    done
    
    while true
    do
        sleep 1
    done
fi
