function [WSC_Max, Ang_Max, DeltaComp_Max] = optimalWSC_ZF(A, E, Rj, hj, dAB, gammaA, gammaJ, channelParam, thetaV, nUAV, typeA )
% Performs exhaustive search of WSC over available actions (opening angles)
% and returns the WSC, angle and Delta distribution corresponding to it for
% the NS precoding case

    WSC_Max = 0;            % Initialize WSC
    Ang_Max = 0;            % Initialize angle
    DeltaComp_Max = 0;      % Initialize Delta distribution
    for ang=thetaV          % Iterate over the available angles (action space search)
        UAVs                    = setNewPos_N(nUAV, ang, hj, Rj, typeA);    % Position UAV in corresponding positions
        [WSC_Dum, DeltaComp]    = computeWSC_ZF_NUAV(A, E, UAVs, dAB, gammaA, gammaJ, channelParam );   % Compute WSC
        
        % If WSC is better than current best WSC, update
        if WSC_Dum > WSC_Max
            WSC_Max = WSC_Dum;
            Ang_Max = ang;
            DeltaComp_Max = DeltaComp;
        end 
    end 
end