#!/bin/bash -l
tglffolder 13;pbsMonitor  -jq parallel01 -jn 1 -cn 13 -wt 00:30:00 -exe  tgyro -e . -n 13
