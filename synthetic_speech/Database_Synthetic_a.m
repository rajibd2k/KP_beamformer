% Voiced - Unvoiced Speech
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all ; clc ; close all ;

mkdir('Database') ;

FS = 8000 ;

dur = 5 ; % s

F0_V_Uv = [ 150 ;  0 ] ;

F_base = [ 500 ; 1500 ; 2500 ; 3500 ; 0 ] ;

F_dev = [ 0 ; 0 ] ;

Bwmid = 250 ;

for F0_num = 1 : length( F0_V_Uv )
    
    F0 = F0_V_Uv(F0_num) ;

    if F0 > 0 
        k = 100 ; % voiced 
        speech_type = 'voiced' ; 
    else
        k = -100 ; % unvoiced
        speech_type = 'unvoiced' ; 
    end

    %---------------------------------------------------------------------------
    % Bandwidths and speech type
    %---------------------------------------------------------------------------

    Bw = Bwmid + k*[-1.5 ; -0.5 ; 0.5 ; 1.5 ; 0] ; % Bandwidth changes from low/high to high/low
    
    F = F_base + F_dev( F0_num ) ;
    [ excitation , glottal , v1 , v2 , v3 , v4 , v5 , vts , lips , speech ] = SFT_Synth_Speech( dur , FS , F0 , F , Bw , speech_type  ) ;
    speech = speech([1000:end]) ; speech = speech - mean(speech) ;  
    
    m = 2 / ( max(speech) - min(speech) ) ;
    c = -( max(speech) + min(speech) ) / ( max(speech) - min(speech) ) ;
    speech = m*speech + c ;
    
    if F0_num == 1
        filename = 'Synth_speech_voiced_SOI' ;
    elseif F0_num == 2
        filename = 'Synth_speech_unvoiced_SOI' ;
    end
    
    save( ['Database/',filename] ,'speech') ;
        
end


intf_type = {'white', 'babble', 'hfchannel'} ;
for idx_intf_type = 1:length(intf_type) 
    
    mkdir(['Database/',intf_type{idx_intf_type}]) ;

    num_samples = length(speech) ;
    signals = zeros( num_samples , 5 ) ;
    
    Intfs = zeros( num_samples , 4 ) ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    addpath('noise') ;
    intf_filename = [intf_type{idx_intf_type}, '.wav'] ;
    orgintf = audioread( intf_filename ) ;
    for idx_intf = 1 : 4
        beg_sample = (idx_intf-1)*num_samples + 1 ;
        end_sample = beg_sample + num_samples - 1 ;
        Intfs(:,idx_intf) = orgintf( beg_sample : end_sample )' ;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    speech = load('Database/Synth_speech_voiced_SOI') ; speech = speech.speech ; % 70 degrees
    signals(:,1) = speech ;
    signals(:,2:5) = Intfs ; % 20 , 30, 130, 160 degrees

    save( ['Database/', intf_type{idx_intf_type}, '/signals_synthetic_1a'] ,'signals') ;

    speech = load('Database/Synth_speech_unvoiced_SOI') ; speech = speech.speech ; % 70 degrees
    signals(:,1) = speech ;
    signals(:,2:5) = Intfs ; % 20 , 30, 130, 160 degrees

    save( ['Database/', intf_type{idx_intf_type}, '/signals_synthetic_2a'] ,'signals') ;

end

M = 2^6 ; % 64 sensors
sensor_noise = randn( 10*8000, M ) ; % noise for 10 s
save(['./Database/sensor_noise'] , 'sensor_noise') ;

clear all ;
