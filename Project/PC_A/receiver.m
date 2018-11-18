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
time_value = 1;
set(audio_recorder,'TimerPeriod',time_value,'TimerFcn',@audioTimerFcn); %use cells to specify function and its arguments

record(audio_recorder); %start recording

% Edit here


%my_recording = getaudiodata(audio_recorder);
end

% CALLBACK FUNCTION
function audioTimerFcn(recObj, event, handles)
    if recObj.UserData.counter < 10
        disp('Slowly!')
        recObj.UserData.counter = recObj.UserData.counter + 1;
    elseif recObj.UserData.counter == 10
        disp('Changing timer value')
        time_value = 0.1;
        set(recObj,'TimerPeriod',time_value); 
        recObj.UserData.counter = recObj.UserData.counter + 1;
    elseif recObj.UserData.counter < 20
        disp('Superfast!')
        recObj.UserData.counter = recObj.UserData.counter + 1;        
    else
        disp('Complete the receiver') 
        stop(recObj);
    end
end
