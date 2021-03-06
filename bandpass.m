function Hd = bandpass(Fs,Fpass1,Fpass2)
%THETA_BANDPASS Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 8.6 and the Signal Processing Toolbox 7.1.
% Generated on: 16-Sep-2015 17:30:26

% Butterworth Bandpass filter designed using FDESIGN.BANDPASS.

mf = min(Fpass1,Fpass2);
if mf < 20 
    fst = 0.5;
elseif mf >= 20 & mf <= 100   
    fst = 5;
else
    fst = 20;
end

% All frequency values are in Hz.
Fstop1 = Fpass1 - fst;           % First Stopband Frequency
Fstop2 = Fpass2 + fst;          % Second Stopband Frequency
Astop1 = 3*fst;%10;             % First Stopband Attenuation (dB)
Apass  = 1;                     % Passband Ripple (dB)
Astop2 = 3*fst; %10;            % Second Stopband Attenuation (dB)
match  = 'stopband';            % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
                      Astop2, Fs);
Hd = design(h, 'butter', 'MatchExactly', match);

