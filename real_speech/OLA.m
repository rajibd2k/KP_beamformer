%%% OverLap and Add (OLA)
%%% signal_frame : signal samples per frame
%%% framesize im ms, frameshift in ms, FS in Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [signal_postprocess , signal] = OLA( signal_frame , framesize_ms, frameshift_ms, FS, total_samples, sample_start)

framesize_samples = framesize_ms * FS / 1000 ;
frameshift_samples = frameshift_ms * FS / 1000 ;

num_frames = size( signal_frame, 2 ) ; 
num_samples = (num_frames - 1)*frameshift_samples + framesize_samples ;

signal = zeros( num_samples , 1 ) ;
for idx_frame = 1 : num_frames
    
    beg_sample = (idx_frame - 1)*frameshift_samples + 1 ;
    end_sample = (idx_frame - 1)*frameshift_samples + framesize_samples ;
    
    signal(beg_sample:end_sample) = signal(beg_sample:end_sample) + signal_frame(:,idx_frame) ;       
end

signal = signal(1:total_samples) ;
% m = 2 / ( max(signal) - min(signal) ) ;
% c = -( max(signal) + min(signal) ) / ( max(signal) - min(signal) ) ;
% signal = m*signal + c ;
signal = signal / max(abs(signal)) ;

signal_postprocess = signal ;
signal_postprocess(1:sample_start) = median(signal_postprocess) ;
% m = 2 / ( max(signal_postprocess) - min(signal_postprocess) ) ;
% c = -( max(signal_postprocess) + min(signal_postprocess) ) / ( max(signal_postprocess) - min(signal_postprocess) ) ;
% signal_postprocess = m*signal_postprocess + c ;
signal_postprocess = signal_postprocess / max(abs(signal_postprocess)) ;

end