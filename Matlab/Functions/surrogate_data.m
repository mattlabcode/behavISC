function [Y] = surrogate_data(X,r)
 
% [Y] = surrogate_data(X,r)
%
% X - time signal (vector)
% r - number of randomizations
%
% Y - time signals with random phases
%
% References:
% Stam CJ (2005) Nonlinear dynamical analysis of EEG and MEG: Review of an
%    emerging field. Clinical Neurophysiology 116:2266-2301.
%
% Description: The script calculates surrogate data from a time series by
% converting it to the frequency domain (fft), randomizing the phase of the
% frequency components, and converting it back to the time domain (ifft).
% --------------------------------------------------------------------------
% B. Herrmann, Email: herrmann.b@gmail.com, 2014-08-04
 
if min(size(X))>1 || length(size(X))>2, error('Error: X needs to be a vector'); end
 
% have a column vector ready
X = X(:);
 
% complex fourier data
C = fft(X,[],1);
 
% amplitude and phase
A = abs(C);
P = angle(C);
 
% number of random points needed
n = round(length(C)/2)-1;
     
% get random phases
rP = rand([n r])*2*pi - pi;
 
% get symmetry in random phases
if ~mod(n,2)
    rP = [repmat(P(1),[1 r]); rP; flipud(rP)*-1];
else
    rP = [repmat(P(1),[1 r]); rP; repmat(P(n+2),[1 r]); flipud(rP)*-1];
end
 
% blow up amplitudes
Am = repmat(A,[1 r]);
 
% calculate new, randomized complex numbers
Cr = Am .* (cos(rP) + 1i.*sin(rP));
 
% convert back to time domain
Y = ifft(Cr,[],1,'symmetric');
