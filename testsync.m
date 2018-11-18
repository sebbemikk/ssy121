clear all
close all
clc
%

%Tx--------------------------------------------------------------------------------------------------------------------------------------
%


fc=1000; % Carrier frequency (dummy)
fsamp = 44000; % sampeling frequency
Tsymb = 0.003; % Symbol period
Tsamp = 1/fsamp; %Sampling period
M = 2; % Bits per symbol
fsymb = 1/Tsymb; % symbol frequency

% Create Pulse shape (RRC: alpha = 0, G = 0.01)
G = 0.003; % Variable for pulse shape
alpha = 0.3; % variable for pulse shape
[tttt,t] = rtrcpuls(alpha, G, fsamp, 4);


databits = [1 1 1 1 1 1];
%databits = randi([0 1],1,10)
prea = [0 1 0 1];
bits = [prea databits]

message = nan(1,1);
i = 1;
while i<= length(bits)/M  % Same as buffer(), puts 1xN bit vector into MxN/M matrix
    a = bits(M*(i-1)+1:M*i);
    message(i,1:M) = a;
    i = i+1;

end

% Map constellation to message
messageDec = bi2de(message,'left-msb') +1; 
const = [1+1i, -1+1i, -1-1i,1-1i];
x = const(messageDec);
    
% Sample up constellation points 
fsfd = (1)*fsamp/fsymb; % Samples per symbol
x_up = upsample(x,fsfd);
   
% Convolve x_up with the pulse shape to create a train of pulses
s = conv(tttt,x_up);

slen = length(s);
tnew=0:Tsamp:(slen-1)*Tsamp;
% Add carrier frequency
carrcos = @(t)sqrt(2)*cos(2*pi*fc*t);
carrsin = @(t)sqrt(2)*sin(2*pi*fc*t);

% Signal to send
ssend = real(s).*carrcos(tnew) + imag(s).*carrsin(tnew);


%
%Rx-------------------------------------------------------------------------------------------------------------------------------------------
%


% Sent signal is received random delay is added (awgn)
pshapeNoDel = downsample(ssend,1);
randDelay = zeros(1,randi([0 100],1));
pshape = awgn([randDelay pshapeNoDel],10);


plen = length(pshape);

preamessage = [];
i = 1;
while i<= length(prea)/M  % Same as buffer(), puts 1xN bit vector into MxN/M matrix
    a = prea(M*(i-1)+1:M*i);
    preamessage(i,1:M) = a;
    i = i+1;

end
 
% Map constellation to message
preamessageDec = bi2de(preamessage,'left-msb') +1; 
preax = const(preamessageDec);

% To make signal of preamble that can be used for xcorr
preambUpSample = upsample(preax,fsamp/fsymb);
preambpshape = conv(preambUpSample, tttt);



tprea=0:Tsamp:(length(preambpshape)-1)*Tsamp;



preacf = real(preambpshape).*carrcos(tprea) + imag(preambpshape).*carrsin(tprea);


autocor = conv(pshape,preacf);
[maxval, ind] = max(autocor)
figure
plot(pshape)
hold on
stem(ind,1)


% Matched filter
mf = conj(fliplr(mf1));
Eh = sum(abs(mf.^2));
pmf = conv(mf,p)./Eh;

% Move signal Ts to left and sample ()
pmf = pmf(length(mf) + Tdel :length(pmf) - length(mf) + 1);

figure
plot(real(pmf),'r')
hold on
plot(imag(pmf),'b')


prealsamp = downsample(real(pmf), fsfd)
pimagsamp = downsample(imag(pmf), fsfd)
stem(1:fsfd:length(pmf),prealsamp)
stem(1:fsfd:length(pmf),pimagsamp)


psamp = prealsamp + 1i.*pimagsamp

scatterplot(psamp)

