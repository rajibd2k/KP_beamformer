%%% OverLap and Add (OLA)
%%% signal_frame : signal samples per frame
%%% framesize im ms, frameshift in ms, FS in Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function signal = OLA( signal_frame , framesize, frameshift, FS )

framesize_samples = framesize * FS / 1000 ;
frameshift_samples = frameshift * FS / 1000 ;

num_frames = size( signal_frame, 2 ) ; 
num_samples = (num_frames - 1)*frameshift_samples + framesize_samples ;

signal = zeros( num_samples , 1 ) ;
for idx_frame = 1 : num_frames
    
    beg_sample = (idx_frame - 1)*frameshift_samples + 1 ;
    end_sample = (idx_frame - 1)*frameshift_samples + framesize_samples ;
    
    signal(beg_sample:end_sample) = signal(beg_sample:end_sample) + signal_frame(:,idx_frame) ;
        
end

m = 2 / ( max(signal) - min(signal) ) ;
c = -( max(signal) + min(signal) ) / ( max(signal) - min(signal) ) ;
signal = m*signal + c ;

end