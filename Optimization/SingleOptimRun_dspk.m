% Currents that modulate spike shape dynamics on multiple timescales induce
% ramping bursts in model respiratory neurons
% Code written by : Dr. Victor Matveev, NJIT

% Main call to optimization method "fminsearch"

clear;
figure;
plotFlag = 1;  % If no plots required, set this to zero

InitParams = [3.5 0.36 260 110 45 2 50 4.8 73 55 60 9 61.15 8.968];

BestParams = fminsearch( @(x) BurstModel_dspk(x, plotFlag), InitParams);

% Compute one last time to make plots:
BurstModel_dspk(BestParams, plotFlag); 

