% calculate the PSD
% Author: Tiberius Georgescu

function [psd,f] = hpsd(A,nos,fs)
N1 = length(A);                     % signal length 
N2 = floor(N1/nos);                 % segment used for psd length

psd0 = 0;
for i = 1:1:nos
    psd0 = psd0 + N2/fs * abs(fft( (A( (i-1)*N2 + 1 : i*N2 )) /N2)).^2;
end
%% Convert two sided PSD to one sided PSD
psd2 = psd0/nos;
psd = psd2( 1 : floor( length(psd2)/2)  +1);
psd(2:end) = 2*psd(2:end); % except for the DC component, the others are doubled
%%
df = fs/N2; % resolution in frequency space
f = (0 : 1 : length(psd)-1)*df;
f = f';
end

% references:
% [1] Michael Cerna and Audrey F Harvey. The fundamentals of fft-based signal analysis
% and measurement. Technical report, Application Note 041, National Instruments,
% 2000.
% [2] Hanspeter Schmid. How to use the fft for signal and noise simulations and
% measurements. University of Applied Sciences Northwestern Switzerland, page 13, 2009.
% [3] Hebestreit, Erik. Thermal properties of levitated nanoparticles. Diss. ETH Zurich, 2017.