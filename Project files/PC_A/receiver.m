function [audio_recorder] = receiver(fc)

fs = 22050; %sampling frequency
audio_recorder = audiorecorder(fs,24,1);% create the recorder

%ADD USER DATA FOR CALLBACK FUNCTION
audio_recorder.UserData.receive_complete = 0;
audio_recorder.UserData.pack  = []; %allocate for data package
audio_recorder.UserData.pwr_spect = []; %allocate for PSD
audio_recorder.UserData.const = []; %allocate for constellation
audio_recorder.UserData.eyed  = []; %allocate for eye diagram
audio_recorder.UserData.counter = 1;
%attach callback function
time_value = 0.01;
set(audio_recorder,'TimerPeriod',time_value,'TimerFcn',@audioTimerFcn); %use cells to specify function and its arguments

record(audio_recorder); %start recording


end

% CALLBACK FUNCTION
function audioTimerFcn(recObj, event, handles)

    if recObj.UserData.counter == 1 %% look if data is coming
        
        disp('looking for data (fast)')
        pshape= getaudiodata(recObj);
        plot(pshape)
        if abs(max(pshape)) > 0.06
            recObj.UserData.counter = recObj.UserData.counter + 1;
             set(recObj,'TimerPeriod',0.5); 
        end
    elseif recObj.UserData.counter < 10 % When |max(sig)| > 0.003 the program slows down and should record each pulse
        disp('receiving data')
        pshape= getaudiodata(recObj);
        plot(pshape)
        recObj.UserData.counter = recObj.UserData.counter + 1;        
    else
        disp('Complete the receiver') 
        stop(recObj);
        %figure
        %X = fftshift(fft(pshape))
        %plot(X)
    end
    
   
end
