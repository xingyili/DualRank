% Usage:
%              >> A_RWRplus(Adj, r_restart, P0, N_max_iter, Eps_min_change, IsNormalized,  NormalizationType);
%              >> [Pt, WAdj ] = A_RWRplus(Adj, r_restart, P0, N_max_iter, Eps_min_change, IsNormalized,  NormalizationType)
% Inputs:
%       AdjSet   = adjacency matrix of networks
%       r_restart   = probability of restart
%       P0   =  initial state of genes
%       N_max_iter   = number of max iteration
%       Eps_min_change   = threshold
%       IsNormalized   =  if the adjacency matrix need to be normalized
%       NormalizationType   = normalization type
% Outputs:
%       Pt  = scores of genes after using network propagation model
%       WAdj  = normalized matrix
function  [Pt, WAdj ]= A_RWRplus(Adj, r_restart, P0, N_max_iter, Eps_min_change, IsNormalized,  NormalizationType)
	%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    if ~exist('N_max_iter','var') || isempty(N_max_iter) || (isnumeric( N_max_iter) && N_max_iter<=1 ) 
        N_max_iter =100; 
    elseif ~isnumeric( N_max_iter)  
        error('N_max_iter should be isnumeric!') ;
    end
    %
    if ~exist('Eps_min_change','var') || isempty(Eps_min_change) 
        Eps_min_change =10^-6; 
    elseif isnumeric( Eps_min_change) && Eps_min_change>=1 
        warning('The Eps_min_change is nomenaning. Reset Eps_min_change to be 10^-6.'); 
        Eps_min_change =10^-6;  
    elseif ~isnumeric( Eps_min_change)  
        error('Eps_min_change should be isnumeric!') ;
    end
    
    if ~exist('IsNormalized','var') || isempty(IsNormalized) 
        IsNormalized = false;  % Adj has been normalized for fast run.   
    end
    
    if ~exist('NormalizationType','var') || isempty(NormalizationType) 
        NormalizationType = 'ProbabilityNormalizationColumn';     
    end        
	% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    
    if IsNormalized 
        WAdj = Adj; 
    else
        % WAdj = getNormalizedMatrix(Adj, 'col', true );
        switch NormalizationType 
            case { 'row' ,'ProbabilityNormalizationRow'} 
                % label propagation with memory    
                WAdj = getNormalizedMatrix(Adj, 'row', true );  
                
            otherwise
                error(['NormalizationType is wrong: ',char( string(NormalizationType) )]); 
        end        
    end
   
    P0 = reshape( P0,[], 1); 
    Pt = P0;
    for T = 1: N_max_iter
        Pt1 = (1-r_restart)*WAdj*Pt + r_restart*P0;
        if all( sum( abs( Pt1-Pt )) < Eps_min_change )
            break;
        end
        Pt = Pt1;
    end
    
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
end
