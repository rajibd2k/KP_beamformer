% Organize TIMIT Database
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all ; clc ; close all ;
mkdir('Database') ;
 
database_folder = '/media/rajib/Windows/Users/rajib/Desktop/MATLAB_Programs/TIMIT/8kHz' ;

dir1 = {'TRAIN';'TEST'} ;

speech_size = cell(1,2) ;
speech = [] ;
sizes = [] ;
for idx1 = 1:length(dir1) 
    name2 = [database_folder , '/' , dir1{idx1}] ;
    dir2 = dir(name2) ;
    dir2 = {dir2.name} ;
    dir2 = dir2(3:end) ;
    
    for idx2 = 1 : length(dir2)
        name3 = [name2 , '/', dir2{idx2}] ;
        dir3 = dir(name3) ;
        dir3 = {dir3.name} ;
        dir3 = dir3(3:end) ;
        
        for idx3 = 1 : length(dir3)
            name4 = [name3 , '/', dir3{idx3}] ;
            dir4 = dir(name4) ;
            dir4 = {dir4.name} ;
            dir4 = dir4(3:end) ;
            
            for idx4 = 1 : length(dir4)
                name5 = [name4 , '/', dir4{idx4}] ;               
                check = strcmp(name5(end-2:end) , 'wav') ;
                check = check + 0 ;
                
                if check == 1
                    speech = [speech ; {audioread(name5)}] ;
                    sizes = [sizes ; length(audioread(name5))] ;
                end

            end %idx4

        end %idx3
         
    end %idx2
 
end %idx1

[ sizes , sizeind ] = sort( sizes, 'ascend' ) ;
speech = speech( sizeind ) ;

sizeind = find( sizes >= 25000 ) ; % ~3s
sizes = sizes(sizeind) ;
speech = speech( sizeind ) ;

speech_size{1} = speech ;
speech_size{2} = sizes ;

save('./Database/speech_size', 'speech_size') ;

speech_size = load('./Database/speech_size') ; 
speech_size = speech_size.speech_size ;
speech = speech_size{1} ;
sizes = speech_size{2} ;

num_speech = length( sizes ) ;
signals_indices = randi([101 num_speech], 4, 100) ;
signals_indices = [1:100 ; signals_indices ] ;

signals_sizes = sizes(signals_indices) ;

intf_type = {'white', 'babble', 'hfchannel'} ;
for idx_intf_type = 1:length(intf_type) 
    
    mkdir(['Database/', intf_type{idx_intf_type}]) ;
    
    for idx_sig = 1 : 20%100 

        signals_size = signals_sizes(1,idx_sig) ;
        signals = zeros( signals_size , 5);

        signal = speech{signals_indices(1,idx_sig)} ;
        signals(:,1) = signal(1:signals_size) ;
        
        Intfs = zeros( signals_size , 4 ) ;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        addpath('noise') ;
        intf_filename = [intf_type{idx_intf_type}, '.wav'] ;
        orgintf = audioread( intf_filename ) ;
        for idx_intf = 1 : 4
            beg_sample = (idx_intf-1)*signals_size + 1 ;
            end_sample = beg_sample + signals_size - 1 ;
            Intfs(:,idx_intf) = orgintf( beg_sample : end_sample )' ;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        signals(:,[2:5]) = Intfs ;

        M = 2^6 ; % 64 sensors
        c = 340 ;
        delta = 10^-2 ;
        FS = 8000 ;
        theta_d = 70 ; theta_u = [20, 30, 130, 160] ; theta = [theta_d theta_u] ;

        signals = SignalsSensorArray(signals, theta, FS, M, delta, c) ;

        sensor_noise = randn( signals_size, M ) ;
        signals(:,6,:) = sensor_noise ;

        save(['./Database/', intf_type{idx_intf_type}, '/signals_', num2str(idx_sig)] , 'signals') ;

    end
end

clear all ; 
%exit ;