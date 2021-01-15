%% Evaluate the Misalignment
%% Target : column vector
%% Outputs : column vectors of outputs

function [ MLA , Mean_Error, Var_Error ] = Misalignment(Outputs , Target)

num_Outputs = size(Outputs,2) ;

Target = Target * ones(1,num_Outputs) ;

Error = Outputs - Target ;

Mean_Error = nanmean( abs(Error) ) ;

Var_Error = var( Error ) ;

MLA = Var_Error ./ var( Target ) ;
MLA = 100 * sqrt(MLA) ;

% MLA = nansum( Error.^2 ) ;
% MLA = MLA ./ nansum( Target.^2 ) ;
% MLA = 100 * sqrt(MLA) ;

end