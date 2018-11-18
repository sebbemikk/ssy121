function transmitter(packet,fc)
    
fsamp = 35000; % sampeling frequency
Tsymb = 0.006; % Symbol period
Tsamp = 1/fsamp; %Sampling period
M = 2; % Bits per symbol(QPSK)
fsymb = 1/Tsymb; % symbol frequency

info_bits=packet';
preamble=[1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,0,0,1,1,0,0,1,1];% bits added to detect targeted signal
bits=[preamble,info_bits];

message=buffer(bits,2)';%transform bits into message
% Map constellation to message
messageDec = bi2de(message,'left-msb') +1;
constellation = [1+1i, -1+1i, -1-1i,1-1i]/sqrt(2);
x = constellation(messageDec);

% Sample up constellation points to ai*d(-kn)
fsfd = (1)*fsamp/fsymb; % Samples per symbol
x_up = upsample(x,fsfd);

% Create Pulse shape
alpha = 0.6; % variable for pulse shape
G = 0.006; % Variable for pulse shape, which determines the symbol time, bandwidth is (1+alpha)/(2*G)
span=4;% which determines number of peaks in the rrc_pulse
[rrc_p, t] = rtrcpuls(alpha,G,fsamp,span); % RRC with alpha and G factors

% Convolve x_up with the pulse shape to create a train of pulses
pulse_train = conv(rrc_p,x_up);

% ================== carrier function================================
carr_cos =@(t) sqrt(2)*cos(2*pi*fc*t);
carr_sin =@(t) sqrt(2)*sin(2*pi*fc*t);
carrier_len = length(pulse_train);
t_new=0:Tsamp:(carrier_len-1)*Tsamp;%define sampling time of carrier

%put baseband signal on carrier
sig_send = real(pulse_train).*carr_cos(t_new) + imag(pulse_train).*carr_sin(t_new);

soundsc(sig_send,fsamp)%send out the signal by soundcard

end