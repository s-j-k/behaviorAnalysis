Fs = 8192; % this is the default sampling frequency in Hz
Ts = 1/Fs; % sampling interval
T=0:Ts:(Fs*Ts);

tone = [1000 4756.82846 5656.854249 6727.171322 ...
    8000 9513.65692 11313.7082 ...
    13454.34264 16000 22627.417];

Y0=zeros(1,Fs*2); % silent interval
Freq=[];
for ii = 1:length(tone)
    Freq(ii,:) = sin(2*pi*tone(ii)*T);
end

Ys=[Freq(1,:) Y0 Freq(2,:) Y0 Freq(3,:) Y0...
    Freq(4,:) Y0 Freq(5,:) Y0 Freq(6,:) Y0...
    Freq(7,:) Y0 Freq(8,:) Y0];
soundsc(Ys, Fs);
