% Finding all the local maxima and minima points
% Input Y is a column vector

function [y_maxima , x_maxima ,y_minima , x_minima ]=extrema(Y)

num_datapoints = size(Y,1) ;

% ---------------------------------------------------------------------------------------------
% calculating slope
 % ---------------------------------------------------------------------------------------------
       
 slope_Y = sign(diff(Y)) ; 
 slope_Y = [0;slope_Y] ; % we consider the slope of the first point to be zero
 
 % ---------------------------------------------------------------------------------------------
 
 % ---------------------------------------------------------------------------------------------
% calculating extrema for case 1
% ----------------------------------------------------------------------------------------------
Y1 = [0;slope_Y] ;
Y2 = [slope_Y;0] ;
diff_slope = Y1 - Y2 ;
diff_slope =  diff_slope(2 : end) ;

x_maxima1 = find(diff_slope == 2) ;
x_minima1 = find(diff_slope == -2) ;

clear Y1 , clear Y2 ; clear diff_slope ;
% --------------------------------------------------------

% ----------------------------------------------------------------------------------------------
% calculating extrema for case 2
% ----------------------------------------------------------------------------------------------
Y1 = [0;slope_Y;0] ;
Y2 = [0;0;slope_Y] ;
Y3 = [slope_Y;0;0] ;

diff_slope1 = Y2 + Y3 ; % to find out where sum of the slopes of the just previous and the just next sample is zero
diff_slope1 =  diff_slope1(2 : (end-1)) ; 

x_extrema_tmp = find(diff_slope1 == 0); 

            % checking if 1st and last point are included in extrema , and removing them
            % --------------------------------------------------------------------------------------------------------
             xfirst = find(x_extrema_tmp == 1) ;
            xlast = find(x_extrema_tmp == num_datapoints) ;

            if length(xlast) ~= 0
                x_extrema_tmp(xlast) = [];
            end
            
            if length(xfirst) ~= 0
                x_extrema_tmp(xfirst) = [];
            end
            
            % --------------------------------------------------------------------------------------------------------
            
y_slope_extrema_tmp = slope_Y(x_extrema_tmp);
x_slope_extrema_tmp = x_extrema_tmp(find(y_slope_extrema_tmp == 0) ) ;

Y_slope_prevsample = slope_Y(x_slope_extrema_tmp - 1) ; % slopes of the just previous sample

Y_slope_indices_maxima = find(Y_slope_prevsample == 1) ;
Y_slope_indices_minima = find(Y_slope_prevsample == -1) ;

x_maxima2 = x_slope_extrema_tmp(Y_slope_indices_maxima) ;
x_minima2 = x_slope_extrema_tmp(Y_slope_indices_minima) ;
            

% -------------------------------------------------------------------------------------------------

x_maxima = sort([x_maxima1;x_maxima2]) ;
x_minima = sort([x_minima1;x_minima2]) ;

y_maxima = Y(x_maxima) ;
y_minima = Y(x_minima) ;

x_maxima = x_maxima - 1 ; % to get the indices starting from 0
x_minima = x_minima - 1 ; % to get the indices starting from 0

clear x_maxima1 , clear x_maxima2 ; clear x_minima1 ; clear x_minima2 ;

end