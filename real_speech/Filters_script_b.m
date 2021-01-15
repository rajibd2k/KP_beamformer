%%% Compute statistics and filters for adaptive beamforming
%%% M = 2^6 ; 
%%% iSIR_dB = 0, iSNR_dB = 10
%%% delta = 10^(-2) ; % 1 cm (inter-sensor-distance)
%%% c = 340 ; % m/s
%%% Ts = 1/8000 ; % s
%%% theta_d = 70 ; theta_u = [20, 30, 130, 160] ; 
%%% iSNR_dB = 10 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all ; close all ; clc ;

addpath(genpath('./pesq-mex-master')) ;

% Datasets
%***********************************************************************************************************************
try
intf_type = {'white', 'babble', 'hfchannel'} ;
for idx_intf_type = 1:length(intf_type) 
        
    for idx_dataset = 1:20 %% 1:100

        data_dir = ['Data/', intf_type{idx_intf_type}] ; %%
        filters_dir = ['Filters_', intf_type{idx_intf_type}, '/', num2str(idx_dataset)] ; %%
        mkdir(filters_dir) ; 

        M = 2^6 ;
        delta = 10^(-2) ; c = 340 ; FS = 8000 ; Ts = 1/FS ;

        framesize_ms = 10 ; % 10 ms
        frameshift_ms = framesize_ms /  2 ; % 5 ms
        FS = 8000 ; Ts = 1/FS ;
        framesize = framesize_ms * FS / 1000 ; 
        frameshift = frameshift_ms * FS / 1000 ; 

        SOI = load([data_dir,'/SOI_' , num2str(idx_dataset)]) ; SOI = SOI.SOI(:,:,1) ; SOI = reshape( SOI, size(SOI,1) , 1 ) ;

        % VAD
        %----------------------------------------------------------------------------
        [ SpeechFrames , SilenceFrames, ~ ] = SFD( SOI, framesize_ms, frameshift_ms, FS , .01 ) ;
        FirstSpeechFrame = SpeechFrames(1) ;
        FirstSpeechSample = (FirstSpeechFrame - 1)*frameshift + 1 ; 


        for iSIR_dB = -10:5:20 

            if sign( iSIR_dB ) == -1
                tmp = 'neg';
            else
                tmp = '' ;
            end
            postname = ['_', num2str(idx_dataset), '_iSNR_dB_10_iSIR_dB_', tmp, num2str(abs(iSIR_dB))] ;
            %----------------------------------------------------------------------------
            Input = load( [data_dir,'/Input', postname] ) ; 
            Input = Input.Input ;

            % Frames / Snapshots 
            %----------------------------------------------------------------------------
            num_samples = size(Input,1) ;
            num_snapshots = ceil( num_samples /  frameshift ) ;
            beg_frames = frameshift*[0:(num_snapshots-1)]' + 1 ;
            end_frames = beg_frames + framesize - 1 ;
            end_frames( find(end_frames > num_samples) ) = num_samples ;

            nfft = 256 ;
            f = [0 : (nfft-1)]' / nfft ;

            Input = reshape( Input, num_samples, M ) ;

            % Compute statistics and filters
            %******************************************************************************************************************************************
            phi = 0:1:180 ; phi_rad = phi/180*pi ;

            z1 = zeros( framesize , num_snapshots ) ; beta_f_1 = zeros( length(phi) , nfft ) ;
            z2 = z1 ; beta_f_2 = beta_f_1;
            z3 = z1 ; beta_f_3 = beta_f_1 ;
            z4 = z1 ; beta_f_4 = beta_f_1 ;
            z5 = z1 ; beta_f_5 = beta_f_1 ;

            phiV = zeros( M, M, length(f) ) ; 

            for idx_snapshots = 1:num_snapshots    

                alpha = 0.8 ; %% speech is dynamic nonstationary signal

                beg_sample = beg_frames(idx_snapshots) ;
                end_sample = end_frames(idx_snapshots) ;

                Y = Input( [beg_sample:end_sample] , : ) ;
                Y = fft( Y , nfft ) ;
                Y = transpose(Y) ;

                for idx_f = 1 : length(f)

                    if idx_snapshots == 1 
                        phiV(:,:,idx_f) = Y(:,idx_f)*Y(:,idx_f)';
                    elseif sum(find( idx_snapshots == SilenceFrames ))
                        phiV(:,:,idx_f) = alpha*phiV(:,:,idx_f) + (1-alpha)*Y(:,idx_f)*Y(:,idx_f)' ; 
                    end

                end

                varname = ['Data/d_M_', num2str(M)] ;
                d = load(varname) ; d = d.d ; %d = d.^Ts ;
                [h] = DS_error( d ) ;
                Z = transpose( nansum( conj(h).*Y ) ) ; z_snapshot = real( ifft( Z ) ) ; 
                z1(:,idx_snapshots) = z_snapshot( 1 : framesize ) ;
                [~, ~, PP_f_norm] = PP(h , delta, c, Ts) ; 
                beta_f_1 = (1-1/idx_snapshots)*beta_f_1  + (1/idx_snapshots)*PP_f_norm ;


                [h] = MVDR_error( phiV , d ) ;
                Z = transpose( nansum( conj(h).*Y ) ) ; z_snapshot = real( ifft( Z ) ) ; 
                z2(:,idx_snapshots) = z_snapshot( 1 : framesize ) ;
                [~, ~, PP_f_norm] = PP(h , delta, c, Ts) ; 
                beta_f_2 = (1-1/idx_snapshots)*beta_f_2  + (1/idx_snapshots)*PP_f_norm ;


                m2 = 3 ; M2 = 2^m2 ; M1 = M / M2 ; n_iter = 5 ;
                varname = ['Data/d_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2) ] ;
                Steering_vectors = load(varname) ; d_1 = Steering_vectors.d_1 ; d_2 = Steering_vectors.d_2 ; %d_1 = d_1.^Ts ; d_2 = d_2.^Ts ;
                [h] = MVDR_Kronecker_error( phiV , d_1 , d_2 , n_iter ) ;
                Z = transpose( nansum( conj(h).*Y ) ) ; z_snapshot = real( ifft( Z ) ) ; 
                z3(:,idx_snapshots) = z_snapshot( 1 : framesize ) ;
                [~, ~, PP_f_norm] = PP(h , delta, c, Ts) ; 
                beta_f_3 = (1-1/idx_snapshots)*beta_f_3  + (1/idx_snapshots)*PP_f_norm ;


                m2 = 3 ; M2 = 2^m2 ; M1 = M / M2 ; n_iter = 1 ;
                varname = ['Data/d_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2) ] ;
                Steering_vectors = load(varname) ; d_1 = Steering_vectors.d_1 ; d_2 = Steering_vectors.d_2 ; %d_1 = d_1.^Ts ; d_2 = d_2.^Ts ;
                [h] = DS_MVDR1_Kronecker_error( phiV , d_1 , d_2 , n_iter ) ;
                Z = transpose( nansum( conj(h).*Y ) ) ; z_snapshot = real( ifft( Z ) ) ; 
                z4(:,idx_snapshots) = z_snapshot( 1 : framesize ) ;
                [~, ~, PP_f_norm] = PP(h , delta, c, Ts) ; 
                beta_f_4 = (1-1/idx_snapshots)*beta_f_4  + (1/idx_snapshots)*PP_f_norm ;

                m2 = 3 ; M2 = 2^m2 ; M1 = M / M2 ; n_iter = 1 ;
                varname = ['Data/d_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2) ] ;
                Steering_vectors = load(varname) ; d_1 = Steering_vectors.d_1 ; d_2 = Steering_vectors.d_2 ; %d_1 = d_1.^Ts ; d_2 = d_2.^Ts ;
                [h] = DS_MVDR2_Kronecker_error( phiV , d_1 , d_2 , n_iter ) ;
                Z = transpose( nansum( conj(h).*Y ) ) ; z_snapshot = real( ifft( Z ) ) ; 
                z5(:,idx_snapshots) = z_snapshot( 1 : framesize ) ;
                [~, ~, PP_f_norm] = PP(h , delta, c, Ts) ; 
                beta_f_5 = (1-1/idx_snapshots)*beta_f_5  + (1/idx_snapshots)*PP_f_norm ;

            end 

            total_samples = num_samples ;  sample_start = FirstSpeechSample ;

            postname = ['_OutputVADbasedEstCurrentDisturbanceStatistics_iSNR_dB_10_iSIR_dB_', tmp, num2str(abs(iSIR_dB))] ;
           %----------------------------------------------------------------------------

            filter = 'DS_F' ;
            [Output_postprocess , Output] = OLA( z1 , framesize_ms, frameshift_ms, FS, total_samples, sample_start) ;
            [ MLA , Mean_Error, Var_Error ] = Misalignment([Output_postprocess , Output] , SOI) ;
            PESQ_NB_BB = [ pesq_mex( SOI, Output, FS, 'both'),  pesq_mex( SOI, Output_postprocess, FS, 'both')]; 
            STOI = [stoi(SOI, Output, FS), stoi(SOI, Output_postprocess, FS)];
            save([filters_dir,'/', filter , postname], 'z1', 'beta_f_1', 'Output', 'Output_postprocess', 'MLA', 'Mean_Error', 'Var_Error', 'PESQ_NB_BB', 'STOI') ; 
%             audiowrite( [filters_dir,'/', filter , postname , '.wav'], Output, FS );
%             audiowrite( [filters_dir,'/', filter , postname1 , '.wav'], Output_postprocess, FS );

            filter = 'MVDR_F' ;
            [Output_postprocess , Output] = OLA( z2 , framesize_ms, frameshift_ms, FS, total_samples, sample_start) ;
            [ MLA , Mean_Error, Var_Error ] = Misalignment([Output_postprocess , Output] , SOI) ;
            PESQ_NB_BB = [ pesq_mex( SOI, Output, FS, 'both'),  pesq_mex( SOI, Output_postprocess, FS, 'both')]; 
            STOI = [stoi(SOI, Output, FS), stoi(SOI, Output_postprocess, FS)];
            save([filters_dir,'/', filter , postname], 'z2', 'beta_f_2', 'Output', 'Output_postprocess', 'MLA', 'Mean_Error', 'Var_Error', 'PESQ_NB_BB', 'STOI') ; 
%             audiowrite( [filters_dir,'/', filter , postname , '.wav'], Output, FS );
%             audiowrite( [filters_dir,'/', filter , postname1 , '.wav'], Output_postprocess, FS );

            filter = 'MVDR_K_F' ;
            [Output_postprocess , Output] = OLA( z3 , framesize_ms, frameshift_ms, FS, total_samples, sample_start) ;
            [ MLA , Mean_Error, Var_Error ] = Misalignment([Output_postprocess , Output] , SOI) ;
            PESQ_NB_BB = [ pesq_mex( SOI, Output, FS, 'both'),  pesq_mex( SOI, Output_postprocess, FS, 'both')]; 
            STOI = [stoi(SOI, Output, FS), stoi(SOI, Output_postprocess, FS)];
            save([filters_dir,'/', filter , postname], 'z3', 'beta_f_3', 'Output', 'Output_postprocess', 'MLA', 'Mean_Error', 'Var_Error', 'PESQ_NB_BB', 'STOI') ; 
%             audiowrite( [filters_dir,'/', filter , postname , '.wav'], Output, FS );
%             audiowrite( [filters_dir,'/', filter , postname1 , '.wav'], Output_postprocess, FS );

            filter = 'DS_MVDR_F' ;
            [Output_postprocess , Output] = OLA( z4 , framesize_ms, frameshift_ms, FS, total_samples, sample_start) ;
            [ MLA , Mean_Error, Var_Error ] = Misalignment([Output_postprocess , Output] , SOI) ;
            PESQ_NB_BB = [ pesq_mex( SOI, Output, FS, 'both'),  pesq_mex( SOI, Output_postprocess, FS, 'both')]; 
            STOI = [stoi(SOI, Output, FS), stoi(SOI, Output_postprocess, FS)];
            save([filters_dir,'/', filter , postname], 'z4', 'beta_f_4', 'Output', 'Output_postprocess', 'MLA', 'Mean_Error', 'Var_Error', 'PESQ_NB_BB', 'STOI') ; 
%             audiowrite( [filters_dir,'/', filter , postname , '.wav'], Output, FS );
%             audiowrite( [filters_dir,'/', filter , postname1 , '.wav'], Output_postprocess, FS );

            filter = 'MVDR_DS_F' ;
            [Output_postprocess , Output] = OLA( z5 , framesize_ms, frameshift_ms, FS, total_samples, sample_start) ;
            [ MLA , Mean_Error, Var_Error ] = Misalignment([Output_postprocess , Output] , SOI) ;
            PESQ_NB_BB = [ pesq_mex( SOI, Output, FS, 'both'),  pesq_mex( SOI, Output_postprocess, FS, 'both')]; 
            STOI = [stoi(SOI, Output, FS), stoi(SOI, Output_postprocess, FS)];
            save([filters_dir,'/', filter , postname], 'z5', 'beta_f_5', 'Output', 'Output_postprocess', 'MLA', 'Mean_Error', 'Var_Error', 'PESQ_NB_BB', 'STOI') ; 
%             audiowrite( [filters_dir,'/', filter , postname , '.wav'], Output, FS );
%             audiowrite( [filters_dir,'/', filter , postname1 , '.wav'], Output_postprocess, FS );

        end

    end

end

exit ;

catch
    
exit ;
    
end
