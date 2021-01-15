%%% Compute actual statistics for adaptive beamforming
%%% M = 2^6 ; 
%%% num_realizations : average statistics must be calculated at these particular realizations
%%% iSIR_dB_range = [-10:5:20]
%%% Statistics for lower values of M can be obtained from subsets of matrices derived for M = 2^6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all ; close all ; clc ;
warning ('off','all') ;

theta_d = 70 ;
theta_u = [20, 30, 130, 160] ; 

delta = 10^(-2) ; % 1 cm (inter-sensor-distance)
c = 340 ; % m/s
Ts = 1/8000 ; % s

iSNR_dB = 10 ;
M = 2^6 ; % num_sensors = 64

num_Intfs = 4 ;
stats_dir = 'Statistics_true_synthetic_1a' ; %%
mkdir(stats_dir) ; 

intf_type = {'white', 'babble', 'hfchannel'} ;
for idx_intf_type = 1:length(intf_type) 

    % These parameters / variables / datasets are constant for all the cases
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    Signals = load(['Data_Snapshots/', intf_type{idx_intf_type}, '/Signals_synthetic_1a']) ; %%
    SOI = Signals.Signals(:,1,:) ;
    Intfs = Signals.Signals(:,[2:(num_Intfs+1)],:) ;
    clear Signals ;

    Sensor_Noise = load('Data_Snapshots/Sensor_Noise') ; 
    Sensor_Noise = Sensor_Noise.Sensor_Noise ; Sensor_Noise = Sensor_Noise(:,[1:M],:) ;

    SOI = reshape( SOI, size(SOI,1), size(SOI,3) ) ;
    %****************************************************************************
    varname = [stats_dir,'/SOI'] ;
    %****************************************************************************
    save(varname, 'SOI') ;

    num_snapshots = size(SOI,2) ;
    snapshots = [1, 10:10:90, 100:50:450, 500:100:num_snapshots, num_snapshots] ;

    nfft = 256 ;
    f = [0 : (nfft-1)]' / nfft ;

    % Variances of SOI and Noise at the first sensor
    var_x = nanmean( var(SOI) ) ;
    var_w = nanmean( nanmean( var( Sensor_Noise ) ) ) ; % noise at all the sensors have identical distribution

    
    % Recalibrate Noise amplitudes based on iSNR
    % Need to be done once, as sensor noise is white only
    %----------------------------------------------------------------------------
    if idx_intf_type == 1
        iSNR = 10^( 0.1 * iSNR_dB ) ;
        tmp_var_w = var_x / iSNR ;
        w_mulfactor = sqrt( tmp_var_w / var_w ) ; 

        %****************************************************************************
        varname = [stats_dir,'/Sensor_Noise_iSNR_dB_10'] ;
        %****************************************************************************
        Sensor_Noise = w_mulfactor * Sensor_Noise ;
        save(varname, 'Sensor_Noise') ;
    end
    

    % Recalibrate interference amplitudes based on iSIR
    %----------------------------------------------------------------------------
    for iSIR_dB = [-10:5:20]

        if sign( iSIR_dB ) == -1
            tmp = 'neg';
        else
            tmp = '' ;
        end

        % Recalibrate interference amplitudes based on iSIR 
        %----------------------------------------------------------------------------
        iSIR = 10^( 0.1 * iSIR_dB ) ;
        tmp_var_u = var_x / iSIR ;
        % Variances of Interferences at the first sensor
        var_u = nanmean( var( Intfs ) , 3  ) ;

        for idx_intf = 1 : num_Intfs
            u_mulfactor = sqrt( tmp_var_u / (num_Intfs * var_u(idx_intf) ) ) ; % assuming the interferences are independent
            Intfs(:,idx_intf,:) = u_mulfactor * Intfs(:,idx_intf,:) ; % hence need to reload data for each iSIR value
        end
        %****************************************************************************
        postname = ['_iSIR_dB_', tmp, num2str(abs(iSIR_dB))] ;
        %****************************************************************************
        mkdir([stats_dir,'/', intf_type{idx_intf_type}]) ;
        save([stats_dir,'/', intf_type{idx_intf_type},'/Intfs', postname], 'Intfs') ; 
        clearvars iSIR tmp_var_u u_mulfactor ;
    end

end

% Sensor positions and Steering vectors
%----------------------------------------------------------------------------
Delta = delta*[0:(M-1)]' ;
varname = [stats_dir,'/Delta_M_', num2str(M) ] ;
save( varname, 'Delta') ;

for m2 = 1:(log2(M)-1)
    
    m1 = log2(M) - m2 ; M1 = 2^m1 ;
    M2 = 2^m2 ;
    
    % actual and virtual steering vectors for the SOI
    [ Steering_vectors , Sensor_positions ] = d1_d2( Delta, M2, theta_d, f, c, Ts ) ;
    d = Steering_vectors.d ; d_1 = Steering_vectors.d_1 ; d_2 = Steering_vectors.d_2 ; 
    varname = [stats_dir,'/d_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2) ] ;
    save( varname, 'd', 'd_1' , 'd_2') ;
    varname = [stats_dir,'/Delta_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2) ] ;
    Delta1 = Sensor_positions.Delta1 ; Delta2 = Sensor_positions.Delta2 ;
    save( varname, 'Delta', 'Delta1', 'Delta2') ;
    
    if m2==1
        varname = [stats_dir,'/d_M_', num2str(M)] ;
        save( varname, 'd') ;
    end
    
    % steering vectors for the interference signals
    if m2==1
        num_Intfs = size(Intfs,2) ;
        du = zeros( M, length(f), num_Intfs ) ;
        for idxintf = 1 : num_Intfs
            [ Steering_vectors , ~ ] = d1_d2( Delta, M2, theta_u(idxintf), f, c, Ts ) ;
            du(:,:,idxintf) = Steering_vectors.d ; 
        end
        varname = [stats_dir,'/du_M_', num2str(M) ] ;
        save( varname, 'du') ;
    end
    
end


% Compute actual statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try   

for idx_intf_type = 1:length(intf_type) 
    
    for iSIR_dB = [-10:5:20]

        if sign( iSIR_dB ) == -1
            tmp = 'neg';
        else
            tmp = '' ;
        end

        %****************************************************************************
        postname = ['_iSIR_dB_', tmp, num2str(abs(iSIR_dB))] ;
        %****************************************************************************
        Intfs = load([stats_dir,'/', intf_type{idx_intf_type},'/Intfs', postname], 'Intfs') ; 
        Intfs = Intfs.Intfs ;

        phiX = zeros( M, M, length(f) ) ;  
        phiW = zeros( M, M, length(f) ) ; 
        phiU = zeros( M, M, length(f) ) ; 

        for idx_snapshots = 1:num_snapshots    

            alpha = 1 - 1/idx_snapshots ;

            fft_SOI = fft(SOI(:,idx_snapshots) , nfft) ;
            fft_Intfs = fft(Intfs(:,:,idx_snapshots) , nfft) ;
            fft_Noise = fft(Sensor_Noise(:,:,idx_snapshots) , nfft) ;

            for idx_f = 1 : length(f)

                    X = d(:,idx_f)*fft_SOI(idx_f) ; 
                    W = fft_Noise(idx_f, :)' ; 

                    U = zeros(M,num_Intfs) ; tmp_phiU = zeros( M, M ) ;
                    for idx_Intfs = 1 : num_Intfs
                         U(:,idx_Intfs) = du(:,idx_f,idx_Intfs) * fft_Intfs(idx_f, idx_Intfs) ;  
                         tmp_phiU = tmp_phiU + U(:,idx_Intfs)*U(:,idx_Intfs)';
                    end

                    if idx_snapshots == 1 
                        phiX(:,:,idx_f) = X*X';
                        phiW(:,:,idx_f) = W*W';
                        phiU(:,:,idx_f) = tmp_phiU ; 
                    else
                        phiX(:,:,idx_f) = alpha*phiX(:,:,idx_f) + (1-alpha)*X*X' ; 
                        phiW(:,:,idx_f) = alpha*phiW(:,:,idx_f) + (1-alpha)*W*W' ;
                        phiU(:,:,idx_f) = alpha*phiU(:,:,idx_f) + (1-alpha)*tmp_phiU ;
                    end

            end

            clear tmp_phiU ;

           if ~sum( find( snapshots == idx_snapshots ) ) 
                continue ;
            end

            postname = ['_M_', num2str(M) , '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10_', 'snapshots_', num2str(idx_snapshots) ] ;
            %----------------------------------------------------- 
            phiV = phiU + phiW ;
            phiY = phiX + phiV ;
            save([stats_dir,'/', intf_type{idx_intf_type},'/phiY', postname], 'phiY') ;         

        end


        if iSIR_dB == -10
            %****************************************************************************
            postname = ['_M_', num2str(M)] ;
            %****************************************************************************   
            save([stats_dir,'/', intf_type{idx_intf_type},'/phiX', postname], 'phiX') ;

            %****************************************************************************
            postname = '' ;
            %****************************************************************************   
            phiX1 = phiX(1,1,:) ; phiX1 = reshape( phiX1, 1, size(phiX1,3) )' ;
            save([stats_dir,'/', intf_type{idx_intf_type},'/phiX1', postname], 'phiX1') ;

            %****************************************************************************
            postname = ['_M_', num2str(M), '_iSNR_dB_10'] ;
            %****************************************************************************   
            save([stats_dir,'/', intf_type{idx_intf_type},'/phiW', postname], 'phiW') ;
        end

        %****************************************************************************
        postname = ['_M_', num2str(M), '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10'] ;
        %****************************************************************************
        save([stats_dir,'/', intf_type{idx_intf_type},'/phiY', postname], 'phiY') ;         
        save([stats_dir,'/', intf_type{idx_intf_type},'/phiV', postname], 'phiV') ;

        phiV1 = phiV(1,1,:) ; phiV1 = reshape( phiV1, 1, size(phiV1,3) )' ;

        iSINR_f = abs( phiX1 ./ phiV1 ) ;
        iSINR_dB_f = 10*log10( iSINR_f ) ;
        iSINR = abs( nansum(phiX1) / nansum(phiV1) ) ;
        iSINR_dB = 10*log10( iSINR ) ;

        %****************************************************************************
        postname = ['_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10'] ;
        %****************************************************************************

        save([stats_dir,'/', intf_type{idx_intf_type},'/phiV1', postname], 'phiV1') ;
        save([stats_dir,'/', intf_type{idx_intf_type},'/iSINR_dB_f', postname], 'iSINR_dB_f') ;
        save([stats_dir,'/', intf_type{idx_intf_type},'/iSINR_dB', postname], 'iSINR_dB') ;      

    end

end
   
warning ('on','all');
exit ;
    
catch
    
warning ('on','all');
exit ;
    
end
