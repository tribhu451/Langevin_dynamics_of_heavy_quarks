#!/bin/bash

nohup ./langevin > run_langevin.log &
JID=$!
wait $JID
echo -e "Langevin run ended at :  "
date

touch "run_finished_flag.dat"



