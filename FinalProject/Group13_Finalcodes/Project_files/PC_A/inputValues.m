classdef inputValues
    properties
        fs
        fc
        fsymb
        pulseShape
        constellation
        pre_carrier
        
    end
    methods
        function inputval = inputValues(fs, fc,fsymb,pulseShape,constellation,pre_carrier)
            inputval.fs = fs;
            inputval.fc = fc;
            inputval.fsymb = fsymb;
            inputval.pulseShape = pulseShape;
            inputval.constellation = constellation;
            inputval.pre_carrier = pre_carrier;
            
        end
    end
end