function sw=spectral_whitening(signal) 
%  Spectral whitening: flattens the spectrum by dividing the foruier 
%                         transform of the signal by the smoothed amplitude
%                         spectrum. 
%  Input: 
%         data = input signal
%
%  Output: 
%       sw = spectrally whitened signal
%
% Written by Karina Løviknes 
% 

L=length(signal); % Length of the signal

fsig=fft(signal,L); % Fourier transform of the signal

% THe amplitude spectrum:   
mag = abs(fsig);

% Smoothening filter:
b=(1/3)*ones(1,3);  % third degree running average filter
magm=filtfilt(b,1,mag);

% A small constant (0.1 % of the amplitude spectrum) is added to avoid 
% dividing  by very small numbers:
prosent=0.1;
wl=prosent*(max(magm))/100;

% The spectral whitening is done by dividing the foruier transform of the 
% signal by the smoothed amplitude spectrum:
swf=fsig./(wl+magm);
sw=ifft(swf,L,'symmetric');

end