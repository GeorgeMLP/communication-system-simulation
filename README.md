# Communication System Simulation

Communication system simulation written in MATLAB. The procedure of our simulation is shown in the following figure.

![Our workflow.](https://github.com/GeorgeMLP/communication-system-simulation/raw/master/Workflow.svg)

First, the sender generates a random binary signal, which is then modulated by BPSK or 64-QAM and oversampled. The resulting signal goes through a low pass filter, and is then added up with additive white Gaussian noise (AWGN).

To restore the noised signal, we first pass the signal through the same low pass filter used previously, and then dawnsample and demodulate the signal, and compare it with the original signal. The error rate is calculated to evaluate the performance of our communication system.

All the ```.m``` files need to be run on MATLAB version R2016b or above, otherwise functions such as ```randi``` or ```rcosdesign``` may report errors.

The content of each file is as follows:

- ```BPSK Using System Functions.m``` adopts BPSK modulation, and makes use of built-in MATLAB functions such as ```rcosdesign```, ```awgn```, etc.
- ```QAM Using System Functions.m``` adopts 64-QAM modulation, and makes use of built-in MATLAB functions such as ```rcosdesign```, ```awgn```, etc.
- ```QAM without System Functions.m``` adopts 64-QAM modulation, and does not use built-in MATLAB functions.
- ```Monte Carlo.m``` exploits the Monte Carlo algorithm to measure the relationship between bit error rate and signal-to-noise ratio.
- ```Comparison between BPSK and 64-QAM.pdf``` compares the error rate of using BPSK and 64-QAM modulation under different signal-to-noise ratios.
- ```Relationship between Error Rate and Signal-to-Noise Ratio.pdf``` plots the relationship between bit error rate and signal-to-noise ratio.
- ```Experiment Report.pdf``` is my experiment report for this assignment, written in Chinese.
- ```Workflow.svg``` shows my workflow.
