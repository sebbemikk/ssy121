function transmitter(packet,fc)


    t1=clock;
    fc=1000; % Carrier frequency (dummy)
    fsamp = 44000; % sampeling frequency
    Tsymb = 0.002; % Symbol period
    Tsamp = 1/fsamp; %Sampling period
    M = 2; % Bits per symbol(QPSK)
    fsymb = 1/Tsymb; % symbol frequency

   


    info_bits = randi([0 1],1,432);
    l_bits=length(info_bits)
    preamble=[1,1,1,1,1,1,1,1,1,1];
    preamble=[preamble,zeros(1,20)];
    bits=[preamble,info_bits];
    message=buffer(bits,2)';
    
    % Map constellation to message
    messageDec = bi2de(message,'left-msb') +1;
    constellation = [1+1i, -1+1i, -1-1i,1-1i]/sqrt(2);
    
    x = constellation(messageDec);

% Sample up constellation points to ai*d(-kn)
    fsfd = (1)*fsamp/fsymb % Samples per symbol
%     fsfd=500;
    x_up = upsample(x,fsfd);
    
 
    
    
    
    subplot(2,1,1)
    stem(0:Tsymb:Tsymb*(length(messageDec)-1),messageDec)
    xlabel('t');
    ylabel('value');
    title('value of bits')
    subplot(2,1,2)
    stem(0:Tsamp:Tsamp*(length(x_up)-1),real(x_up))
    xlabel('t');
    ylabel('value of x_up');
    title('"real" part of upsampled messaged')
    figure   

%     % Create Pulse shape 
%  
    alpha = 0.3; % variable for pulse shape
    G = 0.002; % Variable for pulse shape, which determines the symbol time, bandwidth is (1+alpha)/(2*G)
    span=4;% which determines number of peaks in the rrc_pulse
    [rrc_p, t] = rtrcpuls(alpha,G,fsamp,span); % RRC with alpha and G factors
% % % % % % length of pulse train is 703;    
       
    % Convolve x_up with the pulse shape to create a train of pulses
    pulse_train = conv(rrc_p,x_up);
    


    subplot(2,1,1)
    plot(real(pulse_train)) 
    
    title('real-part of pulse_train(symbol)')
    subplot(2,1,2)
    plot(imag(pulse_train))
    title('imaginary-part of pulse_train(symbol)')
    figure
    
    

    
% ================== carrier function================================
    carr_cos =@(t) sqrt(2)*cos(2*pi*fc*t);
    carr_sin =@(t) sqrt(2)*sin(2*pi*fc*t);
%  ==================================================================    
%  ==================================================================



    carrier_len = length(pulse_train);
    t_new=0:Tsamp:(carrier_len-1)*Tsamp;
   
    sig_send = real(pulse_train).*carr_cos(t_new) + imag(pulse_train).*carr_sin(t_new);
    
    
%   ===========transmitted signal====================================
    plot(t_new,sig_send)
    title('transmitted signal')
    figure
%   =================================================================
%   =================================================================


   
%  ===============frequency of the transmitted signal================  
    SIG_send = abs(fftshift(fft(sig_send)));
    N = length(SIG_send);
    dF = fsamp/N;                      % hertz
    f = -fsamp/2:dF:fsamp/2-dF;
    plot(f,SIG_send)
    title('frequency of transmitted signal')
    figure
%  ==================================================================
%  ==================================================================
    
   
    soundsc(sig_send)
    
    

end