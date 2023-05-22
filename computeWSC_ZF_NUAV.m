function [WSC_Val, DeltaComp] = computeWSC_ZF_NUAV(A, E, UAVs, dAB, gammaA, gammaJ, channelParam )
% Computes the WSC for a given B, UAVs, and considering NS precoding.

    %   A:              Alice's Position        1x3
    %   B:              Bob's Position          1x3
    %   E:              Eve's Positions         nEx3
    %   UAVs:           Position of UAV1        nUAVx3
    %   Rj:             Orbit radius of UAVs    1x1
    %   hj:             UAVs altittude          1x1
    %   dAB:            Distance from A to B	1x3     Can be actual or estimated
    %   gammaA:         Tx SNR at A             1x1
    %   gammaJ:         Tx SNR at UAVs          1x1
    %   channelParam:   Channel and PL parameters
    % [1]   Flores, Alejandro; Moya Osorio, Diana Pamela; Juntti, Markku (2022): A multi-armed bandit 
    %       framework for efficient UAV-based cooperative jamming coverage. TechRxiv. Preprint. 
    %       https://doi.org/10.36227/techrxiv.21564456.v1
    % [2]   X. A. Flores Cabezas, D. P. Moya Osorio and M. Latva-Aho, "Distributed UAV-enabled 
    %       zero-forcing cooperative jamming scheme for safeguarding future wireless networks," 2021 
    %       IEEE 32nd Annual International Symposium on Personal, Indoor and Mobile Radio Communications 
    %       (PIMRC), Helsinki, Finland, 2021, pp. 739-744, doi: 10.1109/PIMRC50174.2021.9569692.
    % [3]   J. P. Vilela, M. Bloch, J. Barros and S. W. McLaughlin, "Wireless Secrecy Regions With 
    %       Friendly Jamming," in IEEE Transactions on Information Forensics and Security, vol. 6, no. 2, 
    %       pp. 256-266, June 2011, doi: 10.1109/TIFS.2011.2111370.
    
    nUAV = size(UAVs,1);                % Obtain nUAV from UAVs (so it's not explicitely passed)
    N = nUAV/2;                         % Number of groups of 2 UAVs
    
    B = A;                              
    B(1) = B(1) + dAB;                  % Locate B with respect to A
    
    % Channel parameters
    phi         = channelParam(1);
    omega       = channelParam(2);
    alpha       = channelParam(3);
    alpha_AG    = channelParam(4);
    ne_LOS      = channelParam(5);
    ne_NLOS     = channelParam(6);
    Rs          = channelParam(7);
    k           = channelParam(8);      % Number of rays of Nakagami channel (m)
    choice      = channelParam(10);     % Choice of which channel to use(0: Rayleigh, 1: Nakagami-m)
    
    % Parameters regarding Eve
    dAE = transpose(sqrt( ( A(:,1) - E(:,1) ).^2 + ( A(:,2) - E(:,2) ).^2 ));   % Distances from A to every E
    OmegaAE = dAE.^alpha;                                                       % Nakagami channel spread

    % UAVs - Es channels
    dJE = sqrt( ( UAVs(:,1) - E(:,1)' ).^2 + ( UAVs(:,2) - E(:,2)' ).^2  + ( UAVs(:,3) - E(:,3)' ).^2); % Distances from UAVs to every E
    Theta_JE = (180/pi) * asin(UAVs(:,3)./dJE);                                                         % Angles from UAVs to every E
    PLOS_JE = 1./(1 + phi * exp( -omega*( Theta_JE - phi ) ) );                                         % Probability of LoS from UAVs to every E
    LJE = PLOS_JE.*(abs(dJE).^alpha_AG)*ne_LOS + (1-PLOS_JE).*(abs(dJE).^alpha_AG)*ne_NLOS;             % Average pathloss from UAVs to every E
    gJE = 1./LJE;               % Deterministic channel gain from UAVs to every E
    hJE = sqrt(gJE);            % Deterministic channel response from UAVs to every E
    
    % Parameters regarding Bob
    OmegaAB = dAB.^alpha;       % Nakagami channel spread
    
    % UAVs - B channels
    dJB = sqrt( ( UAVs(:,1) - B(1) ).^2 + ( UAVs(:,2) - B(2) ).^2  + ( UAVs(:,3) - B(3) ).^2);          % Distances from UAVs to B
    Theta_JB = (180/pi) * asin(UAVs(:,3)./dJB);                                                         % Angles from UAVs to B
    PLOS_JB = 1./(1 + phi * exp( -omega*( Theta_JB - phi ) ) );                                         % Probability of LoS from UAVs to B
    LJB = PLOS_JB.*(abs(dJB).^alpha_AG)*ne_LOS + (1-PLOS_JB).*(abs(dJB).^alpha_AG)*ne_NLOS;             % Average pathloss from UAVs to B
    gJB = 1./LJB;               % Deterministic channel gain from UAVs to B
    hJB = sqrt(gJB);            % Deterministic channel response from UAVs to B

    % Secrecy metrics
    
    % Pair-wise interference gain on E
    g_INT = 0;
    for i=1:N       % Computed per pair of UAVs
       h_INT = hJB(2*i)*hJE(2*i - 1,:) - hJB(2*i - 1)*hJE(2*i,:);           % Pair-wise interference "channel"
       g_INT = g_INT + abs(h_INT).^2;                                       % Pair-wise interference "gain" from (3)
    end
    
    aNJ = gammaA;                           % B parameter for SOP calculation for no-jamming case
    bNJ = gammaA;                           % E parameter for SOP calculation for no-jamming case
    
    aJ = gammaA;                            % B parameter for SOP calculation for jamming case
    bJ = gammaA ./ ( 1 + gammaJ*g_INT );    % E parameter for SOP calculation for jamming case

    % Secrecy calculations
    
    % SOP calculations (for all Es)
    beta    = ((2.^Rs)-1)./aJ;          % Parameter for SOP
	eta     = (2.^Rs).*(bJ./aJ);        % Parameter for SOP
    switch choice       % Choose G2G fading type (0: Rayleigh, 1: Rician)
        case 0      % G2G Rayleigh fading SOP [1]
            SOP_J   = 1 - (exp( -(OmegaAB./aJ ).*(2^Rs - 1) ))*( 1 ./ ( (2^Rs)*(OmegaAB./OmegaAE).*( bJ./aJ ) + 1 ) );
            SOP_NJ  = 1 - (exp( -(OmegaAB./aNJ).*(2^Rs - 1) ))*( 1 ./ ( (2^Rs)*(OmegaAB./OmegaAE).*(bNJ./aNJ) + 1 ) );
        case 1      % G2G Nakagami-m fading SOP [2]
            SOP_J   = SOP_NakagamiM_N(beta,             eta,    1./OmegaAB, 1./OmegaAE, k, k ,1);
            SOP_NJ  = SOP_NakagamiM_N((2.^Rs-1)/gammaA, 2.^Rs,  1./OmegaAB, 1./OmegaAE, k, k ,1);
    end
    
    % Area-based secrecy metrics [1-3]
    DeltaComp = (1-SOP_J)./(1-SOP_NJ);          % Calculation of Delta [1,2]
    coverage    =   sum( DeltaComp(:)>1);       % Jamming coverage [1-3]
    efficiency  =   mean( DeltaComp(:) );       % Jamming efficiency [1-3]
    
    WSC_Val     =   coverage*efficiency;        % Weighted coverage [1,2]

end