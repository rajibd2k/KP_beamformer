%%% Dataset for adaptive beamforming
%%% M = 2^6 ; 
%%% iSIR_dB = [-10:5:20], iSNR_dB = 10
%%% delta = 10^(-2) ; % 1 cm (inter-sensor-distance)
%%% c = 340 ; % m/s
%%% Ts = 1/8000 ; % s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all ; clc ; close all ;
warning ('off','all') ;

theta_d = 70 ;
theta_u = [20, 30, 130, 160] ; 

delta = 10^(-2) ; % 1 cm (inter-sensor-distance)
c = 340 ; % m/s
FS = 8000 ; % Hz
Ts = 1/FS ; % s

iSNR_dB = 10 ;
M = 2^6 ; % num_sensors = 64

num_Intfs = 4 ;
database_dir = 'Database' ; 
data_dir = 'Data' ; 
mkdir(data_dir) ; 


% Datasets
%***********************************************************************************************************************
intf_type = {'white', 'babble', 'hfchannel'} ;
for idx_intf_type = 1:length(intf_type) 
    
    mkdir(['Data/', intf_type{idx_intf_type}]) ;
    
    for idx_dataset = 1:20 %% 1:100

        signals = load([database_dir, '/', intf_type{idx_intf_type}, '/signals_' , num2str(idx_dataset)]) ; 

        SOI = signals.signals(:,1,[1:M]) ; 
        save([data_dir, '/', intf_type{idx_intf_type},'/SOI_' , num2str(idx_dataset)], 'SOI' ) ;
        %audiowrite( [data_dir, '/', intf_type{idx_intf_type},'/original_' , num2str(idx_dataset), '.wav'], SOI(:,:,1), FS );

        Intfs = signals.signals(:,[2:(num_Intfs+1)],[1:M]) ;

        Sensor_Noise = signals.signals(:,end,[1:M]) ;

        clear signals ;

        % Frames / Snapshots 
        %----------------------------------------------------------------------------
        framesize = 10 * FS / 1000 ; % 10 ms
        frameshift = framesize /  2 ; % 5 ms
        num_samples = size(SOI,1) ;
        num_snapshots = ceil( num_samples /  frameshift ) ;
        beg_frames = frameshift*[0:(num_snapshots-1)]' + 1 ;
        end_frames = beg_frames + framesize - 1 ;
        end_frames( find(end_frames > num_samples) ) = num_samples ;

        nfft = 256 ;
        f = [0 : (nfft-1)]' / nfft ;

        % Variance of SOI 
        %----------------------------------------------------------------------------
        var_x = nanmean( var(SOI) , 3 ) ; 

        % Recalibrate signal amplitudes based on iSNR
        %----------------------------------------------------------------------------
        var_w = nanmean( var(Sensor_Noise) , 3 ) ;
        iSNR = 10^( 0.1 * iSNR_dB ) ;
        tmp_var_w = var_x / iSNR ;
        w_mulfactor = sqrt( tmp_var_w / var_w ) ; 
        Sensor_Noise = w_mulfactor * Sensor_Noise ;

        % Recalibrate interference amplitudes based on iSIR
        %----------------------------------------------------------------------------
        for iSIR_dB = [-10:5:20]

            if sign( iSIR_dB ) == -1
                tmp = 'neg';
            else
                tmp = '' ;
            end

            iSIR = 10^( 0.1 * iSIR_dB ) ;
            tmp_var_u = var_x / iSIR ;

            var_u = nanmean( var( Intfs ) , 3  ) ; % Variances of Interferences calculated at each loop

            for idx_intf = 1 : num_Intfs
                u_mulfactor = sqrt( tmp_var_u / (num_Intfs * var_u(idx_intf) ) ) ; % assuming the interferences are independent
                Intfs(:,idx_intf,:) = u_mulfactor * Intfs(:,idx_intf,:) ; % hence need to reload data for each iSIR value
            end

            Intf_SensorNoise = nansum([Intfs , Sensor_Noise],2) ; 
            Input = nansum([SOI , Intfs , Sensor_Noise],2) ; 

            postname = ['_', num2str(idx_dataset), '_iSNR_dB_10_iSIR_dB_', tmp, num2str(abs(iSIR_dB))] ;
            %----------------------------------------------------------------------------
            % save([data_dir, '/', intf_type{idx_intf_type},'/Intf_SensorNoise', postname], 'Intf_SensorNoise') ; 
            save([data_dir, '/', intf_type{idx_intf_type},'/Input', postname], 'Input') ; 
            % audiowrite( [data_dir, '/', intf_type{idx_intf_type},'/Input_first_sensor', postname, '.wav'], Input(:,:,1), FS ) ;

            clearvars iSIR tmp_var_u u_mulfactor ;
        end

        clearvars Intf_SensorNoise Input SOI Intfs Sensor_Noise 

    end

end


% Sensor positions and Steering vectors
%***********************************************************************************************************************

Delta = delta*[0:(M-1)]' ;
varname = [data_dir,'/Delta_M_', num2str(M) ] ;
save( varname, 'Delta') ;

for m2 = 1:(log2(M)-1)

    m1 = log2(M) - m2 ; M1 = 2^m1 ;
    M2 = 2^m2 ;

    % actual and virtual steering vectors for the SOI
    [ Steering_vectors , Sensor_positions ] = d1_d2( Delta, M2, theta_d, f, c, Ts ) ;
    d = Steering_vectors.d ; d_1 = Steering_vectors.d_1 ; d_2 = Steering_vectors.d_2 ; 
    varname = [data_dir,'/d_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2) ] ;
    save( varname, 'd', 'd_1' , 'd_2') ;
    varname = [data_dir,'/Delta_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2) ] ;
    Delta1 = Sensor_positions.Delta1 ; Delta2 = Sensor_positions.Delta2 ;
    save( varname, 'Delta', 'Delta1', 'Delta2') ;

    if m2==1
        varname = [data_dir,'/d_M_', num2str(M)] ;
        save( varname, 'd') ;
    end

    % steering vectors for the interference signals
    if m2==1
        du = zeros( M, length(f), num_Intfs ) ;
        for idx_intf = 1 : num_Intfs
            [ Steering_vectors , ~ ] = d1_d2( Delta, M2, theta_u(idx_intf), f, c, Ts ) ;
            du(:,:,idx_intf) = Steering_vectors.d ; 
        end
        varname = [data_dir,'/du_M_', num2str(M) ] ;
        save( varname, 'du') ;
    end

end

exit ; 