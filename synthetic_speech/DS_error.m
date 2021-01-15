%% Delay-and-Sum (DS) Filter Beamforming 
%% est_d : Steering vector of the SOI at each frequency

function [h] = DS_error( est_d )

h = est_d / size(est_d , 1) ;

end
