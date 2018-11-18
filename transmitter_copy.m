% function transmitter(packet,fc)
    clear all
    close all
    clc

    t1=clock;
    fc=1000; % Carrier frequency (dummy)
    fsamp = 44000; % sampeling frequency
    Tsymb = 0.002; % Symbol period
    Tsamp = 1/fsamp; %Sampling period
    M = 2; % Bits per symbol(QPSK)
    fsymb = 1/Tsymb; % symbol frequency

   

%     info_bits = [1,1,1,0,0,1,0,0,1,1];
%     info_bits = [1 1 1 1 1 1 0 0 1 0 0 1 1 1 0 1 1 1 0 0 1 0 0 1 1 1 0 1 1 1 0 0 1 0 0 1 1 1 1 1 1 1 0 0 1 0 0 1 1 1 0 1 1 1 0 0 1 0 0 1 1 1 0 1 1 1 0 0 1 0 0 1]; % Dummy string of bits
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
%     scatterplot(x); hold on
%     figure()

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
    
   
    
    
    
    
    
%  ========================Reciever Part============================   
%  ========================Reciever Part============================
    delay_data=zeros(1,randi([10,500],1,1));
%     delay_data=zeros(1,260);
    delay_data=awgn(delay_data,10);
    l_delay=length(delay_data)
    sig_receive=[delay_data, sig_send];%%%get a random delayed data
    
    p_message=buffer(preamble,2)'; % transform preamble into message(2 bits per symble)
    mess_Dec = bi2de(p_message,'left-msb') +1;
    constellation = [1+1i, -1+1i, -1-1i,1-1i]/sqrt(2);
    p_mess_up = constellation(mess_Dec);
    
    up_preamble=upsample(p_mess_up,fsamp/fsymb);% upsample the presample-message
    pulse_preamble=conv(up_preamble,rrc_p);
    
    l_pulse_preamble=length(pulse_preamble);
    
    
    %%%%========= only use part of the presample to eliminate =================
    %%%%===============influence caused by bits_message=========================
    l_used=floor(3*l_pulse_preamble/5); %premable samples that are used in synchrolization
    l_left=ceil(2*l_pulse_preamble/5);% premable samples that are not used in synchrolization
    
    pulse_preamble=pulse_preamble(1:l_used);% use part of the preamble
    
    
%   multiply pulse_preamble with carrire so that we can do atuocorrelation
    carr_len = length(pulse_preamble);
    tt_new=0:Tsamp:(carr_len-1)*Tsamp;
    pre_carrier = real(pulse_preamble).*carr_cos(tt_new) + imag(pulse_preamble).*carr_sin(tt_new);
    
    plot(sig_send,'b')
    hold on
    title('waveform of preamble on carrier')
    plot(pre_carrier,'r')
    figure
    
    
    
    corr = conv(sig_receive, fliplr(pre_carrier ));
    Ecorr = sum(pre_carrier.^2);
    [maxvalue,location] = max(abs(corr));
    figure
    plot(corr/Ecorr)
    hold on
    stem(location,maxvalue/Ecorr)
    maxvalue/Ecorr
    
    last_preamble_pos=location; %%%%%  position of the last preamble bit
    fir_preamble=last_preamble_pos-length(pre_carrier )+1%%% the position of the first preamble sample
    
%     plot((abs(corr)))
%     hold on
%     stem(location,maxvalue)
%     title('sychrolization')
%     figure
    
    
% =============================cut off time_delay========================
% =============================cut off time_delay========================

    l_sig_receive=length(sig_receive)
    sig_receive=sig_receive(fir_preamble:l_sig_receive);% cut off the time_delay part
    
% =======================================================================



    t_receive=0:Tsamp:(length(sig_receive)-1)*Tsamp;
    pulse_real=sig_receive.*carr_cos(t_receive);
    pulse_imag=sig_receive.*carr_sin(t_receive);
    
    
    
%   ===============Fourier Transform so as to pass LPF========================
%   ==========================================================================
    fft_real=fft(pulse_real);
    fft_imag=fft(pulse_imag);
%   ==========================================================================
%   ==========================================================================
    
    
%   =======frequency of the received signal multiplied by carrier============= 
%   ==========================================================================
    REAL = abs(fftshift(fft_real));
    IMAG = abs(fftshift(fft_imag));
    N_S = length(REAL);
    dF = fsamp/N_S                      % hertz
    f = -fsamp/2:dF:fsamp/2-dF;
   
    subplot(2,1,1)
    stem(REAL)
    title('frequency of the recieved signal multiplied by cosin-carrier, which is real part')
    subplot(2,1,2)
    plot(f,IMAG)
    title('frequency of the recieved signal multiplied by sin-carrier,which is imaginary part')
    figure
    
% % ==========================================================================
% % ==========================================================================



% 
% % ===============Matched Filter Receiver===================================

rrc_p;
MFR=conj(fliplr(rrc_p));
E_MFR=sum(abs(MFR.^2));

message_real=conv(MFR,pulse_real)./E_MFR;
message_real=message_real(length(MFR):end-length(MFR)+1);

message_imag=conv(MFR,pulse_imag)./E_MFR;
message_imag=message_imag(length(MFR):end-length(MFR)+1);



NN=fsamp/fsymb
re=downsample(message_real,NN);
im=downsample(message_imag,NN);




subplot(2,1,1)
plot(0:Tsamp:Tsamp*(length(message_real)-1),message_real)
hold on
stem(0:Tsamp*NN:(length(re)-1)*Tsamp*NN,re)
title('real part of received message')
subplot(2,1,2)
plot(0:Tsamp:Tsamp*(length(message_imag)-1),message_imag)
hold on
title('imaginary part of received message')
stem(0:Tsamp*NN:(length(im)-1)*Tsamp*NN,im)
% figure

recovered_message=re+1i*im;
l_preamble=length(preamble);
l_rec_message=length(recovered_message);

recovered_message=recovered_message(l_preamble/2+1:l_rec_message);
scatterplot(recovered_message)

% ===========decoding the signal from constellation ==============
% ================================================================
    d_message=[];
    for i=1:length(recovered_message)
        if real(recovered_message(i))>=0&imag(recovered_message(i))>0
            d_message(i)=1;
        elseif real(recovered_message(i))<0&imag(recovered_message(i))>=0
            d_message(i)=2;
        elseif real(recovered_message(i))<=0&imag(recovered_message(i))<0
            d_message(i)=3;
        else
            d_message(i)=4;
        end
    end
    d_message=d_message-1;
    d_bits=de2bi(d_message,'left-msb')';
    d_bits=d_bits(:)';
    length_d_bits=length(d_bits)
    
         
    BER=length(find(info_bits-d_bits)~=0)
    
% ==================================================================
% ==================================================================

    t2=clock;
    run_time=etime(t2,t1)

% % %  % PERFORM PSD CALCULATIONS IN BASE BAND BEFORE MATCHED FILTERING
% % % 
% % %     [pxx, pwr_spect.f] = pwelch(BB_signal,1024,768,1024, fsamp); % note that pwr_spect.f will be normalized frequencies
% % % 
% % %     pwr_spect.f = fftshift(pwr_spect.f); %shift to be centered around zero
% % % 
% % %     pwr_spect.f(1:length(pwr_spect.f)/2) = pwr_spect.f(1:length(pwr_spect.f)/2) - fsamp;
% % % 
% % %     pwr_spect.p = fftshift(10*log10(pxx/max(pxx))); % shift, normalize and convert PSD to dB
% % % 
% % % % end 
    
    
    
    




