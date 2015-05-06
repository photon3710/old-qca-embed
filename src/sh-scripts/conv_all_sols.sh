#!/bin/bash

python conv_sols.py ../sols/bench/dense/lim/4b-acc/
python conv_sols.py ../sols/bench/dense/lim/4b-mux/
python conv_sols.py ../sols/bench/dense/lim/full-add/
python conv_sols.py ../sols/bench/dense/lim/mem-loop/
python conv_sols.py ../sols/bench/dense/lim/sel-osc/
python conv_sols.py ../sols/bench/dense/lim/ser-add/
python conv_sols.py ../sols/bench/dense/lim/sr-flip-flop/
python conv_sols.py ../sols/bench/dense/lim/xor-gate/

python conv_sols.py ../sols/bench/dense/full/4b-acc/
python conv_sols.py ../sols/bench/dense/full/4b-mux/
python conv_sols.py ../sols/bench/dense/full/full-add/
python conv_sols.py ../sols/bench/dense/full/mem-loop/
python conv_sols.py ../sols/bench/dense/full/sel-osc/
python conv_sols.py ../sols/bench/dense/full/ser-add/
python conv_sols.py ../sols/bench/dense/full/sr-flip-flop/
python conv_sols.py ../sols/bench/dense/full/xor-gate/

