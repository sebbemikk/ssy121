clear all
close all
clc
%

%Tx--------------------------------------------------------------------------------------------------------------------------------------
%


fc=1000; % Carrier frequency (dummy)
fsamp = 44000; % sampeling frequency
Tsymb = 0.003; % Symbol period
%Rb = 440; %bit rate
Tsamp = 1/fsamp; %Sampling period
M = 2; % Bits per symbol
fsymb = 1/Tsymb; % symbol frequency

% Create Pulse shape (RRC: alpha = 0, G = 0.01)
G = 0.003; % Variable for pulse shape
%tsy = -Tsymb*4:Tsamp:Tsymb*4; % time vector for making pulse shape
alpha = 0.3; % variable for pulse shape
%rrc = @(t)((sin(pi.*(1-alpha).*t./G)+(4.*alpha.*t./G).*cos(pi.*(1+alpha).*t./G))./(sqrt(G).*(pi.*t./G).*(1-(4.*alpha.*t./G).^2))); % RRC pulse
%tttt=rrc(tsy); % RRC with alpha and G factors

[tttt,t] = rtrcpuls(alpha, G, fsamp, 4);

%plot(t,tttt)
%figure 

%bits = [0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0];
databits = [1 1 1 1 1 1];

%databits = randi([0 1],1,10)
prea = [0 1 0 1];
bits = [prea databits]


%bits = [1 1 1 1 1 1 0 0 1 0 0 1 1 1 0 1 1 1 0 0 1 0 0 1 1 1 0 1 1 1 0 0 1 0 0 1 1 1 1 1 1 1 0 0 1 0 0 1 1 1 0 1 1 1 0 0 1 0 0 1 1 1 0 1 1 1 0 0 1 0 0 1]; % Dummy string of bits
%bits = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];

sizebits = size(bits)
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
    

    
%t=-Tsymb:Tsamp:Tsymb; %??
    
% Sample up constellation points to ai*d(-kn)
fsfd = (1)*fsamp/fsymb; % Samples per symbol
x_up = upsample(x,fsfd);


%stem(real(x_up))
%figure
    
% Convolve x_up with the pulse shape to create a train of pulses
s = conv(tttt,x_up);

%S = abs(fftshift(fft(s)));
%N = length(S);
%dF = fsamp/N;                      % hertz
%f = -fsamp/2:dF:fsamp/2-dF;
%figure
%plot(f,S)

figure() 
plot(real(s),'r')
hold on
plot(imag(s),'b')
slen = length(s);
tnew=0:Tsamp:(slen-1)*Tsamp;

    
carrcos = @(t)sqrt(2)*cos(2*pi*fc*t);
carrsin = @(t)sqrt(2)*sin(2*pi*fc*t);
%figure 
%plot(tnew,carrcos)

ssend = real(s).*carrcos(tnew) + imag(s).*carrsin(tnew);
%figure
%plot(real(ssend))
%subplot(2,1,2)
%plot(imag(ssend))

%Ssend = abs(fftshift(fft(ssend)));
%N = length(Ssend);
%dF = fsamp/N;                      % hertz
%f = -fsamp/2:dF:fsamp/2-dF;
%figure
%plot(f,ssend)


%
%Rx-------------------------------------------------------------------------------------------------------------------------------------------
%


% Sent signal is received and downsampled to receivers samp freq
fs = 44000;
ds = round(fsamp/fs);


pshapeNoDel = downsample(ssend,ds);
randDelay = zeros(1,randi([0 100],1));%10*fsamp/fsymb);
%pshape = pshapeNoDel;
pshape = awgn([randDelay pshapeNoDel],10);
% matched filter ()

plen = length(pshape);
tnew=0:1/fs:(plen-1)*(1/fs);


[mf1,tmf] = rtrcpuls(alpha, G, fs, 4);
mf = conj(fliplr(mf1));
Eh = sum(abs(mf.^2));
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

preambUpSample = upsample(preax,fs);
preambpshape = conv(preambUpSample, mf1);
tprea=0:Tsamp:(length(preambpshape)-1)*Tsamp;


figure
plot(real(preambpshape))
figure
plot(carrcos(tprea))
preacf = real(preambpshape).*carrcos(tprea) + imag(preambpshape).*carrsin(tprea);
figure
plot(preacf)

autocor = conv(pshape,preacf);
[maxval, ind] = max(autocor)

figure
plot(autocor)

figure
plot(pshape)
hold on
stem(ind,maxval)

% FT of received signal
%Pshape = abs(fftshift(fft(pshape)));
%N = length(Pshape);
%dF = fs/N;                      % hertz
%f = -fs/2:dF:fs/2-dF;
%figure
%plot(f,Pshape)

% Multiply Q with cos and I with sin to move freq to 2xfc and baseband

%p = pshape.*carrcosr + 1i.*pshape.*carrsinr;


% FT of signal after cos multiplication
%P = abs(real(fftshift(fft(p))));
%N = length(P);
%dF = fs/N;                      % hertz
%f = -fs/2:dF:fs/2-dF;
%figure
%plot(f,P)

% FT of signal after sin multiplication
%Psin = abs(real(fftshift(fft(psin))));
%N = length(Psin);
%dF = fs/N;                      % hertz
%f = -fs/2:dF:fs/2-dF;
%figure
%plot(f,Psin)

%Filter with mf

pmf = conv(mf,p)./Eh;
figure
plot(real(pmf),'r')
hold on
plot(imag(pmf),'b')

%Tdel = 0;
%preaRx = [0 1 0 1 0 1 0 1 0 1];

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
%%
% THIS CODE BELOW IS TAKEN FROM JOHAN AS A TEST AND SHOULD !!!NOT!!! BE USED
% IN THE PROJECT----------------------------------------------------------
metric = abs(repmat(psamp.',1,4) - repmat(const, length(psamp), 1)).^2; % compute the distance to each possible symbol
[tmp m_hat] = min(metric, [], 2); % find the closest for each received symbol
m_hat = m_hat'-1   % get the index of the symbol in the constellation

%SER = sum(m-1 ~= m_hat) %count symbol errors
m_hat = de2bi(m_hat, 2, 'left-msb')'; %make symbols into bits
a_hat = m_hat(:)' %write as a vector
%BER = sum(bits ~= a_hat) %count of bit errors
