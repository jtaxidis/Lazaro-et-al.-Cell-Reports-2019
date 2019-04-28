function Hd = lowpass(Fs)

Fc = 60;
Nn = 5;
h = fdesign.lowpass('N,F3dB',Nn,Fc,Fs);
Hd = design(h,'butter');

