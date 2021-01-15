%%% Compute filters for adaptive beamforming
%%% M = 2^6 ; 
%%% varying iterations, varying number of snapshots (frames)
%%% iSIR_dB = 0, iSNR_dB = 10
%%% delta = 10^(-2) ; % 1 cm (inter-sensor-distance)
%%% c = 340 ; % m/s
%%% Ts = 1/8000 ; % s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all ; close all ; clc ;
warning ('off','all') ;

theta_d = 70 ;
theta_u = [20, 30, 130, 160] ;
num_Intfs = 4 ;

iSNR_dB = 10 ;
iSIR_dB = 0 ;

if sign( iSIR_dB ) == -1
    tmp = 'neg';
else
    tmp = '' ;
end

M = 2^6 ; % num_sensors = 64
m2_range = [1 : (log2(M) - 1)] ;
n_iter_range = [0, 1 : 2 : 11] ;

stats_dir = ['Statistics_true_synthetic_1a'] ; %%
filters_dir = ['Filters_synthetic_1a'] ; %%
mkdir(filters_dir) ; 

% Dataset
%****************************************************************************
varname = [stats_dir,'/SOI'] ;
%-----------------------------------------------------
SOI = load(varname) ; SOI = SOI.SOI ; 
num_snapshots = size(SOI,2) ; % number of speech frames


intf_type = {'white', 'babble', 'hfchannel'} ;
try 
    
for idx_intf_type = 1:length(intf_type) 

% True Statistics
%****************************************************************************
phiX1 = load([stats_dir,'/', intf_type{idx_intf_type},'/phiX1']) ; phiX1 = phiX1.phiX1 ;

postname = ['_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10'] ;
%-----------------------------------------------------
phiV1 = load([stats_dir,'/', intf_type{idx_intf_type},'/phiV1', postname]) ; phiV1 = phiV1.phiV1 ;
iSINR_dB = load([stats_dir,'/', intf_type{idx_intf_type},'/iSINR_dB', postname]) ; iSINR_dB = iSINR_dB.iSINR_dB ;

postname = ['_M_', num2str(M), '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10'] ;
%-----------------------------------------------------      
phiV = load([stats_dir,'/', intf_type{idx_intf_type},'/phiV', postname]) ; phiV = phiV.phiV ;

    
% Steering vector
%****************************************************************************
varname = [stats_dir,'/d_M_', num2str(M)] ;
d = load(varname) ; d = d.d ;

varname = [stats_dir,'/du_M_', num2str(M)] ;
du = load(varname) ; du = du.du ;


% Compute statistics and filters
%******************************************************************************************************************************************
mkdir([filters_dir,'/', intf_type{idx_intf_type}]) ;

postname = ['_M_', num2str(M) , '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10' ] ;
%----------------------------------------------------- 
% DS_F
[h] = DS_error( d ) ;
Perf = Metrics(h, d([1:M],:), phiX1, phiV1, phiV([1:M],[1:M],:)) ;
save(['./', filters_dir,'/', intf_type{idx_intf_type}, '/DS_F', postname], 'h', 'Perf') ; 


% Obtain statistics and Compute filters
%*******************************************************************************************
for idx_snapshots = [1, 10:10:90, 100:50:450, 500:100:num_snapshots, num_snapshots]   

    postname = ['_M_', num2str(M) , '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10_', 'snapshots_', num2str(idx_snapshots) ] ;
    %----------------------------------------------------- 

    phiY = load([stats_dir,'/', intf_type{idx_intf_type},'/phiY',postname]) ; phiY = phiY.phiY ;

    % MVDR_F
    [h] = MVDR_error( phiY , d ) ;
    Perf = Metrics(h, d([1:M],:), phiX1, phiV1, phiV([1:M],[1:M],:)) ; 
    save(['./', filters_dir,'/', intf_type{idx_intf_type}, '/MVDR_F', postname], 'h', 'Perf') ; 

    for idx_m2 = 1 : length(m2_range) 

        m2 = m2_range( idx_m2 ) ; M2 = 2^m2 ; M1 = M / M2 ;
        varname = [stats_dir,'/d_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2) ] ;
        Steering_vectors = load(varname) ; d_1 = Steering_vectors.d_1 ; d_2 = Steering_vectors.d_2 ;

        for idx_iter = 1 : length(n_iter_range)

            n_iter = n_iter_range( idx_iter ) ;

            postname = ['_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2), '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10_', 'snapshots_', num2str(idx_snapshots), '_iterations_', num2str(n_iter)] ;
            %-----------------------------------------------------

            % MVDR_K_F
            [h] = MVDR_Kronecker_error( phiY , d_1 , d_2 , n_iter ) ;
            Perf = Metrics(h, d([1:M],:), phiX1, phiV1, phiV([1:M],[1:M],:)) ; 
            save(['./', filters_dir,'/', intf_type{idx_intf_type}, '/MVDR_K_F', postname], 'h', 'Perf') ; 

            % DS_MVDR_F
            [h] = DS_MVDR1_Kronecker_error( phiY , d_1 , d_2 , n_iter ) ;
            Perf = Metrics(h, d([1:M],:), phiX1, phiV1, phiV([1:M],[1:M],:)) ; 
            save(['./', filters_dir,'/', intf_type{idx_intf_type}, '/DS_MVDR_F', postname], 'h', 'Perf') ; 

            % MVDR_DS_F
            [h] = DS_MVDR2_Kronecker_error( phiY , d_1 , d_2 , n_iter ) ;
            Perf = Metrics(h, d([1:M],:), phiX1, phiV1, phiV([1:M],[1:M],:)) ; 
            save(['./', filters_dir,'/', intf_type{idx_intf_type}, '/MVDR_DS_F', postname], 'h', 'Perf') ; 

        end

    end

end
    

end
    
    warning ('on','all');
    exit ;

catch

    warning ('on','all');
    exit ;

end