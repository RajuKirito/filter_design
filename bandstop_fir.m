f_samp = 260e3;

%Band Edge speifications
fp1 = 53.5e3;
fs1 = 57.5e3;
fs2 = 77.5e3;
fp2 = 82.5e3;

%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

Wn = [(fs1+fp1)/2 (fs2+fp2)/2]*2/f_samp;        %average value of the two paramters
N_min = ceil((A-7.95) / (2.285*0.096664418));       %empirical formula for N_min

%Window length for Kaiser Window
n=N_min + 12;

%Ideal bandstop impulse response of length "n"

bs_ideal =  ideal_lp(pi,n) - ideal_lp(0.78*pi,n) + ideal_lp(0.5646154*pi,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandStop = bs_ideal .* kaiser_win;
fvtool(FIR_BandStop);         %frequency response

%magnitude response
[H,f] = freqz(FIR_BandStop,1,1024, f_samp);
plot(f,abs(H))
grid