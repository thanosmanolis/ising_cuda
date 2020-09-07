# Bash Scripts

## Description
Scripts made in order to test the executables with various values of n and k, and store the time elapsed of all the executions in a text file, for further observations.

## How to run

1. Copy all these scripts to the main repo folder.
2. Execution: ``./submit_test.sh arg``

The argument has to be one of the 4 implementations (seq / v1 / v2 / v3).

By default, the script will execute the requested implementation with the following specs:

* n: 1000 - 10000 (with a step of 1000)
* k: 50, 100

If you want to change the specs, you can change the variables min_n, max_n, min_k, max_k e.t.c. inside the bash script.

The time elapsed of each execution will be written in a produced text file.