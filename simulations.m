%% MAB formulation for WSC optimization
clc, clear, close all
verStr = '_v2.mat';      % Name of saved data

nMC = 1e3;          % Number of Monte Carlo iterations
typeA = 2 ;         % Angular distribution | 1: symmetric angles, 2: semicircle angles 

% Parameters **************************************************************

% Environment parameters (Urban)
phi = 9.61;         % Environmental constant
omega = 0.16 ;      % Environmental constant
alpha = 0.3;        % Ground Path Loss exponent
alpha_AG = 0.3;     % Air-to-Ground Path Loss Exponent
ne_LOS = 1.0;       % Air-to-Ground LOS attenuation
ne_NLOS = 20;       % Air-to-Ground NLOS attenuation

% Channel parameters
m = 3;              % Number of parallel channels for Nakagami channel
sigma = 1/sqrt(2);  % Noise std dev of each component (real and imag) of every parallel channel
choice = 1;         % 0: Rayleigh channel, 1: Nakagami-m channel
Rs = 1;
channelParam =  [   phi,...             % Store channel parameters in array for function compacity
                    omega,...
                    alpha,...
                    alpha_AG,...
                    ne_LOS,...
                    ne_NLOS,...
                    Rs,...
                    m,...
                    sigma,...
                    choice];

% Nodes *******************************************************************

% Alice
A = [0,0,0];        %   Position of Alice (zero point)
gammaA = 100;       %   Alice Tx SNR

% Bob
sigAB = 1;          %   Unreliability of B's position

% Eves
nR = 50;            %   Number of radial points
nTheta = 180;       %   Number of angular points
nE = nR*nTheta;     %   Total number of Eves

rLow = 0.1;         %   Lowest radius of Eve
rHigh = 50;         %   Highest radius of Eve (Ra)
thetaLow = 0;       %   Lowest angle of Eve
thetaHigh = 2*pi;   %   Highest angle of Eve

rangeR = linspace(rLow,rHigh,nR);                       % Points in Radial dimension
rangeTheta = linspace(thetaLow,thetaHigh,nTheta);       % Points in Angular dimension

[rt, thetat] = meshgrid(rangeR, rangeTheta);            % Radial-angular mesh

E = [rt(:).*cos(thetat(:)), rt(:).*sin(thetat(:)) , zeros(nR*nTheta,1)];    % Eves' position (rectangle coordinates)

% surf(rt.*cos(thetat),rt.*sin(thetat),zeros(nTheta,nR))    % Check form



% UAVs
nUAV = 6;                           %   Number of simultaneous UAVs
gammaJ = gammaA/nUAV;               %   UAVs Jamming SNR
            
hj = 10;                            % UAVs common fixed altitude
Rj = 30;                            % UAVs common fixed orbit around A



% MAB formulation **********************************************************
nAng = 18;                                              % Angle discretization level (opening angle)  -> Number of Angle Actions
if typeA==1 
    angleUAV    = linspace(0,2*pi/(nUAV-1),nAng);       % Possible angle actions (opening angles) for symmetric openings
elseif typeA==2
    angleUAV    = linspace(0,pi/(nUAV-1),nAng);         % Possible angle actions (opening angles) for semicircle openings
end

nLoops = nAng*4;            % Number of loops for action choosing
initWSC = 0;                % Optimistic initial action values
c = 0.3;                    % Exploration parameter for UCB
alpha = 0.1;                % Step size (0: uniform average)


% Auxiliary Variables *****************************************************

% Timing variables
dt = 0;
alphat = 0.2;

% Result-storing variables
WSC_RL      = zeros(nMC,nLoops);        % Store WSC values obtained by MAB formulation
WSC_Max_ZF  = zeros(nMC,nLoops);        % Store maximum WSC values obtained by MAB formulation

Ang_RL_V        = zeros(nMC,nLoops);    % Store angle obtained by MAB formulation
Ang_Max_ZF_V    = zeros(nMC,nLoops);    % Store maximum angle obtained by MAB formulation

dists_AB = zeros(nMC,nLoops);           % Store positions of B

% /*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*

for iMC =1:nMC
    tic
    % Bob's movement
    xo = 0;                     
    xf = 0;
    while abs(xo-xf)<=nLoops*1e-2       % Do not allow movements that are too short (such that each step is more than 0.01)
        xo = rHigh*rand();              % Initial position
        xf = rHigh*rand();              % Final position
    end
    dAB_R_V = linspace(xo,xf,nLoops);	% Total trajectory for all methods (for same MC iteration)
    
    
    % Exhaustive Search results
    for i=1:nLoops
        dAB_R = dAB_R_V(i);         % Bob's new position
        
        % Obtain the maximum WSC for NS precoding and no precoding case
        % for every position of B
        [WSC_Max_Val, Ang_Max_Val, ~]       = optimalWSC_ZF(A, E, Rj, hj, dAB_R, gammaA, gammaJ, channelParam, angleUAV, nUAV, typeA );
        
        WSC_Max_ZF(iMC,i)            = WSC_Max_Val;
        Ang_Max_ZF_V(iMC,i)          = Ang_Max_Val;
    end
    
    % Initialization for RL process
    WSCEst_Angle    = initWSC*ones(1,nAng);                 %   Action value estimation vector for angle actions
    WSCN_Angle      = zeros(1,nAng);                        %   Vector to store angle action ocurrences

    % Results for MAB over iterations and movement of B
    for i=1:nLoops
        dAB_R = dAB_R_V(i);
        
        % dAB estimation and parameter computation
        dAB = normrnd(dAB_R,sigAB);         %   Estimate of the position of B

        %   RL iteration - NS
        [WSCEst_Angle, WSCN_Angle]      = computeRL_UCB(WSCEst_Angle, hj, Rj, WSCN_Angle, angleUAV, ...
                                            A, E, dAB, gammaA, gammaJ, c, channelParam, i, alpha, nUAV,1, typeA);

        %   True WSC values obtained calculations

        % NS
        [~, Ang_RL_Ind]                 = max(WSCEst_Angle);
        Ang_RL                          = angleUAV(Ang_RL_Ind);
        UAVs                            = setNewPos_N(nUAV, Ang_RL, hj, Rj, typeA);
        WSC_RL(iMC,i)          = computeWSC_ZF_NUAV(A, E, UAVs, dAB_R, gammaA, gammaJ, channelParam );

        
        Ang_RL_V(iMC,i)        = Ang_RL;
        
        dists_AB(iMC,i)        = dAB_R;
    end

    % Timing ----------------------------------------------------------
    t1 = toc;                                                        %|
    if iMC == 1                                                      %|
        dt = t1;                                                     %|
    else                                                             %|
        dt = dt + alphat*(t1 - dt);                                  %|
    end                                                              %|
                                                                     %|
    TT = dt*nMC;                                                     %|
    TTF = TT - (iMC-1)*dt;                                           %|
                                                                     %|
    TTF_S = rem(TTF,60);                                             %|
    TTF_M = rem(fix(TTF/60),60);                                     %|
    TTF_H = fix(fix(TTF/60)/60);                                     %|
    % -----------------------------------------------------------------

    fprintf('MC Loop: %i / %i\t\t   Time per Loop: %.2f s\t\t TTF: %i H %i M %.1f S \n',iMC,nMC,t1,TTF_H,TTF_M,TTF_S);    
end
save(['output_nMC-',num2str(nMC), '_nUAV-',num2str(nUAV), '_hj-', num2str(hj), '_Rj-',num2str(Rj), verStr])


%% Grapher
clc, clear, close all
verStr = 'output_nMC-1000_nUAV-6_hj-30_Rj-10_v2.mat';

load(verStr)

size(mean(WSC_RL,1))

figure
hold on
plot(mean(WSC_RL,1)/nE)
plot(mean(WSC_Max_ZF,1)/nE)

grid minor
ylabel("Normalized WSC")
xlabel("Iteration")
legend("MAB", "ES")

figure
hold on
plot(mean(Ang_RL_V,1)*180/pi)
plot(mean(Ang_Max_ZF_V,1)*180/pi)

grid minor
ylabel("Angular separation [degrees]")
xlabel("Iteration")
legend("MAB", "ES")
%%

%%
%% MAB formulation for WSC optimization (multiple variable testing)
clc, clear, close all

load('benchmark.mat', 'B_vec', 'rJ_vec')               % External data for benchmark computation

hj_V        = [30, 50, 100];        % Iterate variable 1
sigmaAB_V   = [0.1, 1, 10];         % Iterate variable 2


vecMain = hj_V;                     % Iterate variable 1 ambiguation (for generalization)
vecVals = sigmaAB_V;                % Iterate variable 1 ambiguation (for generalization)
nMain = length(vecMain);            % Length of iterate variable 1
nVals = length(vecVals);            % Length of iterate variable 1

strMain = 'hj';
strVals = 'sigmaAB';
texMain = '$h_{\mathrm{J}}$';
texVals = '$\sigma_{\mathrm{AB}}$';
verStr = 'data_output_v0.mat';


nMC = 1e3;          % Number of Monte Carlo iterations
typeA = 2 ;         % Angular distribution | 1: symmetric angles, 2: semicircle angles 

% Parameters **************************************************************

% Environment parameters (Urban)
phi = 9.61;         % Environmental constant
omega = 0.16 ;      % Environmental constant
alpha = 0.3;        % Ground Path Loss exponent
alpha_AG = 0.3;     % Air-to-Ground Path Loss Exponent
ne_LOS = 1.0;       % Air-to-Ground LOS attenuation
ne_NLOS = 20;       % Air-to-Ground NLOS attenuation

% Channel parameters
m = 3;              % Number of parallel channels for Nakagami channel
sigma = 1/sqrt(2);  % Noise std dev of each component (real and imag) of every parallel channel
choice = 1;         % 0: Rayleigh channel, 1: Nakagami-m channel
Rs = 1;
channelParam =  [   phi,...             % Store channel parameters in array for function compacity
                    omega,...
                    alpha,...
                    alpha_AG,...
                    ne_LOS,...
                    ne_NLOS,...
                    Rs,...
                    m,...
                    sigma,...
                    choice];

% Nodes *******************************************************************

% Alice
A = [0,0,0];        %   Position of Alice (zero point)
gammaA = 100;       %   Alice Tx SNR

% Bob
% sigAB = 1;          %   Unreliability of B's position (Commented because it is an iterate variable)

% Eves
nR = 50;            %   Number of radial points
nTheta = 180;       %   Number of angular points
nE = nR*nTheta;     %   Total number of Eves

rLow = 0.1;         %   Lowest radius of Eve
rHigh = 50;         %   Highest radius of Eve (Ra)
thetaLow = 0;       %   Lowest angle of Eve
thetaHigh = 2*pi;   %   Highest angle of Eve

rangeR = linspace(rLow,rHigh,nR);                       % Points in Radial dimension
rangeTheta = linspace(thetaLow,thetaHigh,nTheta);       % Points in Angular dimension

[rt, thetat] = meshgrid(rangeR, rangeTheta);            % Radial-angular mesh

E = [rt(:).*cos(thetat(:)), rt(:).*sin(thetat(:)) , zeros(nR*nTheta,1)];    % Eves' position (rectangle coordinates)

% surf(rt.*cos(thetat),rt.*sin(thetat),zeros(nTheta,nR))    % Check form

% UAVs
nUAV = 2;                           %   Number of simultaneous UAVs
gammaJ = gammaA/nUAV;               %   UAVs Jamming SNR

% hj = 70;                          % UAVs common fixed altitude (commented because it is an iterate variable)
Rj = 30;                            % UAVs common fixed orbit around A



% MAB formulation **********************************************************

nLoops = 20;                % Number of loops for action choosing
initWSC = 0;                % Optimistic initial action values
c = 0.3;                    % Exploration parameter for UCB
alpha = 0.1;                % Step size (0: uniform average)

nAng = 10;                                              % Angle discretization level (opening angle)  -> Number of Angle Actions
if typeA==1 
    angleUAV    = linspace(0,2*pi/(nUAV-1),nAng);       % Possible angle actions (opening angles) for symmetric openings
elseif typeA==2
    angleUAV    = linspace(0,pi/(nUAV-1),nAng);         % Possible angle actions (opening angles) for semicircle openings
end


% Auxiliary Variables *****************************************************

% Timing variables
dt = 0;
alphat = 0.2;

% Result-storing variables
WSC_RL      = zeros(nMC,nLoops,nMain,nVals);        % Store WSC values obtained by MAB formulation
WSC_NOP     = zeros(nMC,nLoops,nMain,nVals);        % Store WSC values obtained by 
WSC_GD      = zeros(nMC,nLoops,nMain,nVals);        % Store WSC values obtained by PGD approach
WSC_Max_ZF  = zeros(nMC,nLoops,nMain,nVals);        % Store maximum WSC values obtained by MAB formulation
WSC_Max_NOP = zeros(nMC,nLoops,nMain,nVals);        % Store WSC values obtained by 

WSC_BM      = zeros(nMC,nLoops,nMain,nVals);        % Store WSC values obtained by benchmark

Ang_RL_V        = zeros(nMC,nLoops,nMain,nVals);    % Store angle obtained by MAB formulation
Ang_NOP_V       = zeros(nMC,nLoops,nMain,nVals);    % Store angle obtained by 
Ang_GD_V        = zeros(nMC,nLoops,nMain,nVals);    % Store angle obtained by PGD formulation
Ang_Max_ZF_V    = zeros(nMC,nLoops,nMain,nVals);    % Store maximum angle obtained by MAB formulation
Ang_Max_NOP_V   = zeros(nMC,nLoops,nMain,nVals);    % Store maximum angle obtained by 

dists_AB = zeros(nMC,nLoops,nMain,nVals);           % Store positions of B

% /*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*

for iMC =1:nMC
    
    % Bob's movement
    xo = 0;                     
    xf = 0;
    while abs(xo-xf)<=nLoops*1e-2       % Do not allow movements that are too short (such that each step is more than 0.01)
        xo = rHigh*rand();              % Initial position
        xf = rHigh*rand();              % Final position
    end
    dAB_R_V = linspace(xo,xf,nLoops);	% Total trajectory for all methods (for same MC iteration)
    
    for iVm = 1:nMain                   % Iteration over main iterate
        hj = vecMain(iVm);              % Iterate variable
        tic
        
        
        % Exhaustive Search results
        for i=1:nLoops
            dAB_R = dAB_R_V(i);         % Bob's new position
            
            % Obtain the maximum WSC for NS precoding and no precoding case
            % for every position of B
            [WSC_Max_Val, Ang_Max_Val, ~]       = optimalWSC_ZF(A, E, Rj, hj, dAB_R, gammaA, gammaJ, channelParam, angleUAV, nUAV, typeA );
            [WSC_Max_ValNOP, Ang_Max_ValNOP, ~] = optimalWSC_NOP(A, E, Rj, hj, dAB_R, gammaA, gammaJ, channelParam, angleUAV, nUAV, typeA );

            WSC_Max_ZF(iMC,i,iVm, :)            = WSC_Max_Val;
            Ang_Max_ZF_V(iMC,i,iVm, :)          = Ang_Max_Val;

            WSC_Max_NOP(iMC,i,iVm, :)           = WSC_Max_ValNOP;
            Ang_Max_NOP_V(iMC,i,iVm, :)         = Ang_Max_ValNOP;
        end
        
        for iVa = 1:nVals                   % Iteration over secondary iterate
            sigAB = vecVals(iVa);           % Iterate variable
            
            % Initialization for RL process
            WSCEst_Angle    = initWSC*ones(1,nAng);                 %   Action value estimation vector for angle actions
            WSCN_Angle      = zeros(1,nAng);                        %   Vector to store angle action ocurrences

            WSCEst_Ang_NOP    = initWSC*ones(1,nAng);               %   Action value estimation vector for angle actions
            WSCN_Ang_NOP      = zeros(1,nAng);                      %   Vector to store angle action ocurrences

            WSC_Ang         = angleUAV(fix(length(angleUAV)/2));    %   Angle initialization
            WSC_UAV = 0;                                            %   WSC initial value

            % Initialization for GD process
            Ang_Step = WSC_Ang;                             %   The same angle as the RL approach, to better compare the two approaches

            % Results for MAB over iterations and movement of B
            for i=1:nLoops
                dAB_R = dAB_R_V(i);
                
                % dAB estimation and parameter computation
                dAB = normrnd(dAB_R,sigAB);         %   Estimate of the position of B

                %   RL iteration - NS
                [WSCEst_Angle, WSCN_Angle]      = computeRL_UCB(WSCEst_Angle, hj, Rj, WSCN_Angle, angleUAV, ...
                                                    A, E, dAB, gammaA, gammaJ, c, channelParam, i, alpha, nUAV,1, typeA);

                %   RL iteration - NOP
                [WSCEst_Ang_NOP, WSCN_Ang_NOP]  = computeRL_UCB(WSCEst_Ang_NOP, hj, Rj, WSCN_Ang_NOP, angleUAV, ...
                                                    A, E, dAB, gammaA, gammaJ, c, channelParam, i, alpha, nUAV,3, typeA);
                
                %   GD iteration
                Ang_Step                        = computeGD(A, E, Rj, hj, dAB, gammaA, gammaJ, channelParam, alpha, Ang_Step, i, nUAV, typeA );
                
                %   Benchmark iteration
                UAV_BM                          = computeBM(dAB, hj, B_vec, rJ_vec );

                %   True WSC values obtained calculations

                % NS
                [~, Ang_RL_Ind]                 = max(WSCEst_Angle);
                Ang_RL                          = angleUAV(Ang_RL_Ind);
                UAVs                            = setNewPos_N(nUAV, Ang_RL, hj, Rj, typeA);
                WSC_RL(iMC,i,iVm, iVa)          = computeWSC_ZF_NUAV(A, E, UAVs, dAB_R, gammaA, gammaJ, channelParam );

                % NOP
                [~, Ang_NOP_Ind]                = max(WSCEst_Ang_NOP);
                Ang_NOP                         = angleUAV(Ang_NOP_Ind);
                UAVs                            = setNewPos_N(nUAV, Ang_NOP, hj, Rj, typeA);
                WSC_NOP(iMC,i,iVm, iVa)         = computeWSC_NOP_NUAV(A, E, UAVs, dAB_R, gammaA, gammaJ, channelParam );

                % GD
                UAVs                            = setNewPos_N(nUAV, Ang_Step, hj, Rj, typeA);
                WSC_GD(iMC,i,iVm, iVa)          = computeWSC_ZF_NUAV(A, E, UAVs, dAB_R, gammaA, gammaJ, channelParam );
                
                % Benchmark
                WSC_BM(iMC,i,iVm, iVa)          = computeWSC_1UAV(A, E, UAV_BM, dAB_R, gammaA, gammaJ, channelParam );

                Ang_RL_V(iMC,i,iVm, iVa)        = Ang_RL;
                Ang_NOP_V(iMC,i,iVm, iVa)       = Ang_NOP;
                Ang_GD_V(iMC,i,iVm, iVa)      = Ang_Step;

                dists_AB(iMC,i,iVm, iVa)        = dAB_R;
            end
        end
        
        % Timing ----------------------------------------------------------
        t1 = toc;                                                        %|
        if iMC*iVm == 1                                                  %|
            dt = t1;                                                     %|
        else                                                             %|
            dt = dt + alphat*(t1 - dt);                                  %|
        end                                                              %|
                                                                         %|
        TT = dt*nMC*nMain;                                               %|
        TTF = TT - ( (iMC-1)*nMain + iVm )*dt;                           %|
                                                                         %|
        TTF_S = rem(TTF,60);                                             %|
        TTF_M = rem(fix(TTF/60),60);                                     %|
        TTF_H = fix(fix(TTF/60)/60);                                     %|
        % -----------------------------------------------------------------

        fprintf('MC Loop: %i / %i\t\t %s: %.3f \t\t %s: %.3f \t\t  Time per Loop: %.2f s\t\t TTF: %i H %i M %.1f S \n',iMC,nMC,strMain,vecMain(iVm),strVals,vecVals(iVa),t1,TTF_H,TTF_M,TTF_S);
    end
    
end
save(['data-',strMain, '-',strVals, '-AngT-', num2str(typeA),'-nMC', num2str(nMC) , '-', verStr ])
