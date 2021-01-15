%%% Compute filters for adaptive beamforming
%%% M = 2^6 ; 
%%% iSIR_dB_range = [-10:5:20], iSNR_dB = 10
%%% delta = 10^(-2) ; % 1 cm (inter-sensor-distance)
%%% c = 340 ; % m/s
%%% Ts = 1/8000 ; % s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all ; close all ; clc ;
warning ('off','all') ;

theta_d = 70 ;
theta_u = [20, 30, 130, 160] ;
num_Intfs = 4 ;

delta = 10^(-2) ; % 1 cm (inter-sensor-distance)
c = 340 ; % m/s
Ts = 1/8000 ; % s

iSNR_dB = 10 ;
iSIR_dB_range = [-10:5:20] ;

M = 2^6 ; % num_sensors = 64

stats_dir = ['Statistics_true_synthetic_1a'] ; %%
filters_dir = ['Filters_synthetic_1a'] ; %%
mkdir(filters_dir) ; 

% Dataset
%****************************************************************************
varname = [stats_dir,'/SOI'] ;
%-----------------------------------------------------
SOI = load(varname) ; SOI = SOI.SOI ; 
num_snapshots = size(SOI,2) ; % number of speech frames
idx_snapshots = 100 ; num_snapshots ; %% [10, 100, num_snapshots]
    
% Steering vector
%****************************************************************************
varname = [stats_dir,'/d_M_', num2str(M)] ;
d = load(varname) ; d = d.d ;

varname = [stats_dir,'/du_M_', num2str(M)] ;
du = load(varname) ; du = du.du ;


intf_type = {'white', 'babble', 'hfchannel'} ;
% Compute statistics and filters
%******************************************************************************************************************************************
try
    
for idx_intf_type = 1:length(intf_type)     
       
    % Obtain statistics and Compute filters
    %*******************************************************************************************
    for idx_iSIR_dB = 1:length(iSIR_dB_range)   
        
        iSIR_dB = iSIR_dB_range(idx_iSIR_dB) ;

        if sign( iSIR_dB ) == -1
            tmp = 'neg';
        else
            tmp = '' ;
        end      

        % True Statistics
        %****************************************************************************
        phiX1 = load([stats_dir,'/', intf_type{idx_intf_type},'/phiX1']) ; phiX1 = phiX1.phiX1 ;

        postname = ['_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10'] ;
        %-----------------------------------------------------
        phiV1 = load([stats_dir,'/', intf_type{idx_intf_type},'/phiV1', postname]) ; phiV1 = phiV1.phiV1 ;

        postname = ['_M_', num2str(M), '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10'] ;
        %-----------------------------------------------------      
        phiV = load([stats_dir,'/', intf_type{idx_intf_type},'/phiV', postname]) ; phiV = phiV.phiV ;
              
        
        
        
        postname = ['_M_', num2str(M) , '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10' ] ;
        %----------------------------------------------------- 
        % DS_F
        [h] = DS_error( d ) ;
        Perf = Metrics(h, d([1:M],:), phiX1, phiV1, phiV([1:M],[1:M],:)) ;
        save(['./', filters_dir,'/', intf_type{idx_intf_type}, '/DS_F', postname], 'h', 'Perf') ; 
        
        
        
        
        postname = ['_M_', num2str(M) , '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10_', 'snapshots_', num2str(idx_snapshots) ] ;
        %----------------------------------------------------- 
        phiY = load([stats_dir,'/', intf_type{idx_intf_type},'/phiY',postname]) ; phiY = phiY.phiY ; 
        
        
        postname = ['_M_', num2str(M) , '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10_', 'snapshots_', num2str(idx_snapshots) ] ;
        %----------------------------------------------------- 
        % MVDR_F
        [h] = MVDR_error( phiY , d ) ;
        Perf = Metrics(h, d([1:M],:), phiX1, phiV1, phiV([1:M],[1:M],:)) ; 
        save(['./', filters_dir,'/', intf_type{idx_intf_type}, '/MVDR_F', postname], 'h', 'Perf') ; 

        
        m2 = 3 ; M2 = 2^m2 ; M1 = M / M2 ; n_iter = 5 ;
        varname = [stats_dir,'/d_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2) ] ;
        Steering_vectors = load(varname) ; d_1 = Steering_vectors.d_1 ; d_2 = Steering_vectors.d_2 ;
        postname = ['_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2), '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10_', 'snapshots_', num2str(idx_snapshots), '_iterations_', num2str(n_iter)] ;
        %-----------------------------------------------------
        % MVDR_K_F
        [h] = MVDR_Kronecker_error( phiY , d_1 , d_2 , n_iter ) ;
        Perf = Metrics(h, d([1:M],:), phiX1, phiV1, phiV([1:M],[1:M],:)) ; 
        save(['./', filters_dir,'/', intf_type{idx_intf_type}, '/MVDR_K_F', postname], 'h', 'Perf') ; 
        

        m2 = 3 ; M2 = 2^m2 ; M1 = M / M2 ; n_iter = 1 ;
        varname = [stats_dir,'/d_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2) ] ;
        Steering_vectors = load(varname) ; d_1 = Steering_vectors.d_1 ; d_2 = Steering_vectors.d_2 ;
        postname = ['_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2), '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10_', 'snapshots_', num2str(idx_snapshots), '_iterations_', num2str(n_iter)] ;
        %-----------------------------------------------------
        % DS_MVDR_F
        [h] = DS_MVDR1_Kronecker_error( phiY , d_1 , d_2 , n_iter ) ;
        Perf = Metrics(h, d([1:M],:), phiX1, phiV1, phiV([1:M],[1:M],:)) ; 
        save(['./', filters_dir,'/', intf_type{idx_intf_type}, '/DS_MVDR_F', postname], 'h', 'Perf') ; 
        

        m2 = 3 ; M2 = 2^m2 ; M1 = M / M2 ; n_iter = 1 ;
        varname = [stats_dir,'/d_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2) ] ;
        Steering_vectors = load(varname) ; d_1 = Steering_vectors.d_1 ; d_2 = Steering_vectors.d_2 ;
        postname = ['_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2), '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10_', 'snapshots_', num2str(idx_snapshots), '_iterations_', num2str(n_iter)] ;
        %-----------------------------------------------------
        % MVDR_DS_F
        [h] = DS_MVDR2_Kronecker_error( phiY , d_1 , d_2 , n_iter ) ;
        Perf = Metrics(h, d([1:M],:), phiX1, phiV1, phiV([1:M],[1:M],:)) ; 
        save(['./', filters_dir,'/', intf_type{idx_intf_type}, '/MVDR_DS_F', postname], 'h', 'Perf') ; 
        
    end
    
end
    
    warning ('on','all');
    exit ;

catch

    warning ('on','all');
    exit ;

end