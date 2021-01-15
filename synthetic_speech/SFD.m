%%% Speech Frames Detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ SpeechFrames , SilenceFrames ] = SFD( SOI )

energy_frame = mean( SOI.^2 ) ;
average_energy = mean( energy_frame ) ;

SpeechFrames = find( energy_frame > 0.01*average_energy ) ;

SilenceFrames = find( energy_frame <= 0.01*average_energy ) ;

end