clc
cla 
clear all
close all 
fsamp = 44000; %[samp/sec]
Rb = 440; % [bits/s]
Tsamp = 1/fsamp; % [sec/samp]
M = 2; %number of bits per symbol
fsymb=Rb/2; %[s/symbol]
Tsymb = 1/fsymb; % 






%Bits -> messages 
bits = [0,0,1,1,0,1,1,0,0,0,0,0,1,1,1,1]; %sting of bits to be transmitted
message = nan(1,1) % vector of bits
i=1;
while i <= length(bits)/M
   a = bits(M*(i-1)+1:M*i);
   message(i,1:M) = a;
   
   i=i+1; 
end
message

messageDec=bi2de(message,'left-msb')+1;

% Message -> Symbol
const = [1+1i,-1+1i,-1-1i,1-1i]/sqrt(2); %constellation (sqrt(2) for QAM)
x = const(messageDec); % takes index of value in in message and assigns const point

fsfd= fsamp/fsymb % [samples/symbol]
x_up = upsample(x,fsfd) % upsampled symbols

%plot(x_up)


%Symbol -> Signal

pulse = ones(1,fsfd) % signal shape

s = conv(pulse,x_up) % concovolve symbols with pulse shape

plot(real(s))
hold on 
plot (imag(s),'r')




