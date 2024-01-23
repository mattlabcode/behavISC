function plot_fft(sig, sf)
% Function plots the power spectrum for a signal.
% INPUT:
%   sig - time series (vector of observations for one variable)
%   sf  - data sampling rate
% EXAMPLE:
%{
    SR  = 150;          % sampling frequency
    t   = 0:1/fs:2*pi;  % time vector 
    f   = 5;            % signal frequency
	sig = sin(2*pi*t*f);
    plot(t, sig)
    plot_fft(sig, SR)
%}
% Author: Matthew Bain

% set amount of zero padding (length of fft)
nfft = 1024;

% compute dft
z = fft(sig, nfft);

% get params
N     = length(z);          % number of samples
f     = (0 : N - 1)*(sf/N); % frequency range
power = abs(z).^2/N;        % power of the DFT

% plot power spectrum as function of frequency
plot(f, power, 'color', 'black')

% fft is symmetric so throw away second half
xlim([0 sf/2])    