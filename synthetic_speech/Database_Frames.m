% Speech 10 ms segments (snapshots) for speech sentences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all ; clc ; close all ;
database = 'Database' ;
data_snapshots = 'Data_Snapshots';
mkdir(data_snapshots) ;

FS = 8000 ;
framesize = FS * 10 / 1000 ; %10 ms
frameshift = framesize/2 ; %5 ms

intf_type = {'white', 'babble', 'hfchannel'} ;
for idx_intf_type = 1:length(intf_type) 
    
    database_dir = dir( ['Database/',intf_type{idx_intf_type}] ) ;
    database_dir = {database_dir.name} ;
    database_dir = database_dir(3:end) ;

    for idx_signals = 1 : length(database_dir)

        name = database_dir{idx_signals} ;
        if ~strcmp(name(1:7),'signals') 
            continue ;
        end

        signals = load(['Database/',intf_type{idx_intf_type},'/',name]) ;
        signals = signals.signals ;

        num_samples = size( signals , 1 ) ;
        num_signals = size( signals , 2 ) ;

        beg_frames = [1:frameshift:num_samples]' ;
        end_frames = beg_frames + framesize - 1 ;
        last_frames = find(end_frames > num_samples) ;
        end_frames(last_frames) = num_samples ;

        num_frames = length(beg_frames) ;
        Signals = zeros(framesize,num_signals,num_frames) ;

        for idx_frames = 1 : num_frames

            samples = [beg_frames(idx_frames) : end_frames(idx_frames)]';
            length_samples = length(samples) ;

            Signals([1:length_samples],:,idx_frames) = signals(samples,:) ;

        end
        
        mkdir( [data_snapshots,'/', intf_type{idx_intf_type}] ) ;
        save([data_snapshots,'/', intf_type{idx_intf_type}, '/S',name(2:end)] , 'Signals') ;

    end

end


database_dir = dir(database) ;
database_dir = {database_dir.name} ;
database_dir = database_dir(3:end) ;
for idx_signals = 1 : length(database_dir)
    
    name = database_dir{idx_signals} ;
    if ~strcmp(name(1:(end-4)),'sensor_noise')
        continue ;
    end
    
    signals = load([database,'/',name]) ;
    signals = signals.sensor_noise ;
    
    num_samples = size( signals , 1 ) ;
    num_signals = size( signals , 2 ) ;
    
    beg_frames = [1:frameshift:num_samples]' ;
    end_frames = beg_frames + framesize - 1 ;
    last_frames = find(end_frames > num_samples) ;
    end_frames(last_frames) = num_samples ;
    
    num_frames = length(beg_frames) ;
    Sensor_Noise = zeros(framesize,num_signals,num_frames) ;

    for idx_frames = 1 : num_frames

        samples = [beg_frames(idx_frames) : end_frames(idx_frames)]';
        length_samples = length(samples) ;

        Sensor_Noise([1:length_samples],:,idx_frames) = signals(samples,:) ;

    end

    save([data_snapshots,'/Sensor_Noise'] , 'Sensor_Noise') ;
    
end