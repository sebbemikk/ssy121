function [audio_recorder] = receiver(fc)

fs = 35000; %sampling frequency
Tsymb = 0.006; % Symbol period
Ts = 1/fs; %Sampling period
M = 2; % Bits per symbol(QPSK)
fsymb = 1/Tsymb; % symbol frequency

alpha = 0.6; % variable for pulse shape
G = 0.006; % Variable for pulse shape, which determines the symbol time, bandwidth is (1+alpha)/(2*G)
span=4;% which determines number of peaks in the rrc_pulse
[rrc_p, t] = rtrcpuls(alpha,G,fs,span); % RRC with alpha and G factors
    
preamble=[1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,0,0,1,1,0,0,1,1];% bits added before packet, so as to detect target siganl
p_message=buffer(preamble,2)'; % transform preamble into message(2 bits per symble)

% Map constellation to message
mess_Dec = bi2de(p_message,'left-msb') +1;
constellation = [1+1i, -1+1i, -1-1i,1-1i]/sqrt(2);
p_mess_up = constellation(mess_Dec);

up_preamble=upsample(p_mess_up,round(fs/fsymb));% upsample the presample-message
pulse_preamble=conv(up_preamble,rrc_p);

l_pulse_preamble=length(pulse_preamble);
l_used=floor(3*l_pulse_preamble/5); %premable samples that are used in synchrolization
l_left=ceil(2*l_pulse_preamble/5);% premable samples that are not used in synchrolization   
pulse_preamble=pulse_preamble(1:l_used);% use part of the preamble, so that the transient of bits behand preamble won't influence convolution

carr_len = length(pulse_preamble);
tt_new=0:Ts:(carr_len-1)*Ts;% define sampling time of carrier
carr_cos =@(t) sqrt(2)*cos(2*pi*fc*t);
carr_sin =@(t) sqrt(2)*sin(2*pi*fc*t);
%%%% put preamble-pulse_train on carrier( passband signal)    
pre_carrier = real(pulse_preamble).*carr_cos(tt_new) + imag(pulse_preamble).*carr_sin(tt_new);

inputval = inputValues(fs,fc, fsymb, rrc_p, constellation, pre_carrier);% store variables in function so that it can be used in audioTimerFcn(recObj, event, ~)


audio_recorder = audiorecorder(fs,24,1);% create the recorder
%ADD USER DATA FOR CALLBACK FUNCTION
audio_recorder.UserData.receive_complete = 0;
audio_recorder.UserData.pack  = []; %allocate for data package
audio_recorder.UserData.pwr_spect = []; %allocate for PSD
audio_recorder.UserData.const = []; %allocate for constellation
audio_recorder.UserData.eyed  = []; %allocate for eye diagram
audio_recorder.UserData.inputval= inputval; % Input vars for the callback function

%attach callback function
time_value = 0.1;
set(audio_recorder,'TimerPeriod',time_value,'TimerFcn',@audioTimerFcn); %use cells to specify function and its arguments
record(audio_recorder); %start recording

end

% CALLBACK FUNCTION
function audioTimerFcn(recObj, event, ~)
pshape= (getaudiodata(recObj)).'; %store recieved audio signal in "pshape" 

Ts = 1/recObj.UserData.inputval.fs; %Sampling period
fc=recObj.UserData.inputval.fc; % Carrier frequency 

% Looking for a convolution that passes threshhold --> data being received
Ecorr = sum(recObj.UserData.inputval.pre_carrier.^2);
corr = conv(pshape, fliplr(recObj.UserData.inputval.pre_carrier))/Ecorr;
[maxval,index] = max(corr);
dataSig = [];

if maxval/Ecorr > 0.05  % If data is received, wait for package and save   
    disp('found data :)')
    pause(0.5)% make sure we recieve the whole targeted signal
    
    pshape = getaudiodata(recObj);%store recieved audio signal in "pshape" 
    corr = conv(pshape, fliplr(recObj.UserData.inputval.pre_carrier))/Ecorr;
    [maxval,index] = max(corr); % find the preamble
    
    startData = index -length(recObj.UserData.inputval.pre_carrier)+1;% start point of targeted signal
    endData = startData+(216+13)*(recObj.UserData.inputval.fs/recObj.UserData.inputval.fsymb)+ length(recObj.UserData.inputval.pulseShape)-1;% end point of targeted signal
    dataSig = pshape(startData:endData);% store the targetd signal
    
    % Passband --> baseband
    tt_new=0:Ts:(length(dataSig)-1)*Ts;%define sampling time
    carr_cos =@(t) sqrt(2)*cos(2*pi*fc*t);
    carr_sin =@(t) sqrt(2)*sin(2*pi*fc*t);
    data_bb = carr_cos(tt_new).*dataSig.' + 1i*carr_sin(tt_new).*dataSig.';% Passband --> baseband
    
    % Matched filter
    Errc = sum(recObj.UserData.inputval.pulseShape.^2);
    data_mf = conv(data_bb, conj(fliplr(recObj.UserData.inputval.pulseShape)))/Errc;
    data_mf = data_mf(length(recObj.UserData.inputval.pulseShape):length(data_mf));% signal processed by match filter, but still have tail trancient
    
    fsfd = round(recObj.UserData.inputval.fs/recObj.UserData.inputval.fsymb);
   
    %Symbol --> samples
    data_samp = downsample(data_mf,fsfd);
    firstPar = data_samp(1);
    data_samp = data_samp(14:216+13);% take the targeted message
    
    % Fix attenuation and angle offset
    angOffset = angle(firstPar) - angle(1 -1i);
    ampOffset = abs(firstPar);
    data_samp = (exp(-1i*angOffset)/ampOffset).*data_samp;
    
    % Samples --> bits
    recovered_message = data_samp;
    d_message=[];
    for ii=1:length(recovered_message)
        if real(recovered_message(ii))>=0&imag(recovered_message(ii))>0
            d_message(ii)=1;
        elseif real(recovered_message(ii))<0&imag(recovered_message(ii))>=0
            d_message(ii)=2;
        elseif real(recovered_message(ii))<=0&imag(recovered_message(ii))<0
            d_message(ii)=3;
        else
            d_message(ii)=4;
        end
    end
    d_message=d_message-1;
    d_bits=de2bi(d_message,'left-msb')';
    d_bits=d_bits(:)';
    
    % Step 1: save the estimated bits
    recObj.UserData.pack = d_bits;
    
    % Step 2: save the sampled symbols
    recObj.UserData.const = data_samp;
    
    % Step 3: provide the matched filter output for the eye diagram
    recObj.UserData.eyed.fsfd = fsfd;
    recObj.UserData.eyed.r = data_mf(1:end-9*fsfd+1);% remove the transcient caused by match filter   
    
    % Step 4: Compute the PSD and save it. Note that it has to be computed on
    % the BASE BAND signal BEFORE matched filtering
    [pxx, f] = pwelch(data_mf,1024,768,1024, recObj.UserData.inputval.fs); % note that pwr_spect.f will be normalized frequencies
    f = fftshift(f); %shift to be centered around fs
    f(1:length(f)/2) = f(1:length(f)/2) - recObj.UserData.inputval.fs; % center to be around zero
    p = fftshift(10*log10(pxx/max(pxx))); % shift, normalize and convert PSD to dB
    recObj.UserData.pwr_spect.f = f;
    recObj.UserData.pwr_spect.p = p;
      
    recObj.UserData.receive_complete = 1;
    stop(recObj)       
end
   
end
