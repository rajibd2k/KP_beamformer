%%% Speech Frames Detection
%%% Signal = signals(:,1,1) ; framesize_ms = 10 ; frameshift_ms = 5 ; FS = 8000 ; threshold = .01 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ SpeechFrames , SilenceFrames, SpeechSilenceFrames ] = SFD( Signal, framesize_ms, frameshift_ms, FS , threshold )

framesize = framesize_ms * FS / 1000 ; 
frameshift = frameshift_ms * FS / 1000 ; 
num_samples = size(Signal,1) ;
num_snapshots = ceil( num_samples /  frameshift ) ;

beg_frames = frameshift*[0:(num_snapshots-1)]' + 1 ;
end_frames = beg_frames + framesize - 1 ;
end_frames( find(end_frames > num_samples) ) = num_samples ;

energy_frames = zeros( num_snapshots , 1 ) ;
for idx_snapshots = 1:num_snapshots
    SignalFrame = Signal( beg_frames(idx_snapshots) : end_frames(idx_snapshots) ) ;
    energy_frames(idx_snapshots) = mean(SignalFrame.^2) ;
end

energy_frames = energy_frames - min(energy_frames) ;
energy_frames = energy_frames / max(energy_frames) ;

SpeechSilenceFrames = energy_frames > threshold  ;
SpeechSilenceFrames = SpeechSilenceFrames + 0 ;
SpeechSilenceFrames = smooth(SpeechSilenceFrames, 5, 'moving') ;
SpeechSilenceFrames = SpeechSilenceFrames > 0 ;
SpeechSilenceFrames = SpeechSilenceFrames + 0 ;

SpeechFrames = find( SpeechSilenceFrames > 0 ) ;
SilenceFrames = find( SpeechSilenceFrames == 0 ) ;

end