% A_RWR_M() -  Network propagation
%
% Usage:
%              >> A_RWR_M(AdjSet, P0, restart, pro_jump, UsedGlobalVar, options, NormalizationType);
%              >> [Pset, TableScores ] = ...
%                    A_RWR_M(AdjSet, P0, restart, pro_jump, UsedGlobalVar, options, NormalizationType);
% Inputs:
%       AdjSet   = adjacency matrix of networks
%       P0   =  initial state of genes
%       restart   = the parameter that quantifies the relative importance of the initial state
%       pro_jump   =  the parameter that denotes the relative importance of the states of the neighbors in the layer or the states of nodes in other layers
%       UsedGlobalVar   = adjacency matrix of networks
%       options   =  []
%       NormalizationType   = normalization type
% Outputs:
%       TableScores   = Scores of genes
function [TableScores] = A_RWR_M(AdjSet, P0, restart, pro_jump, UsedGlobalVar, options, NormalizationType)
    global   Global_Var_RWRM 
    if ~exist('UsedGlobalVar','var') || isempty (UsedGlobalVar)
        UsedGlobalVar = false; 
    end   
    tau_vec =[];
    isdebug = false; 
    if exist('options','var') 
        if isfield(options, 'tau_vec')
            tau_vec = options.tau_vec ;   
        end
        if isfield(options, 'isdebug')
            isdebug =options.isdebug ; 
        end
    end   
    if not( UsedGlobalVar )  || isempty( Global_Var_RWRM ) 
       InterLayerMappingType ='FullMapping'; 
        [A_nLxnL,N_node, L_net] = getMultiplexMatrixFromAdjSet(AdjSet, pro_jump, NormalizationType, InterLayerMappingType ) ;   
        Global_Var_RWRM.A_nLxnL_Normalized = A_nLxnL;  
        Global_Var_RWRM.N_node  = N_node   ;
        Global_Var_RWRM.L_net   = L_net   ;
        Global_Var_RWRM.NxL_net = N_node*L_net   ;
        Global_Var_RWRM.pro_jump = pro_jump  ;
        IsNormalized = true; 
        if isdebug; disp('Set Global_Var_RWRM'); end   
    else 
        IsNormalized = true; 
        N_node  = Global_Var_RWRM.N_node;  
        L_net   = Global_Var_RWRM.L_net;  
        pro_jump   = Global_Var_RWRM.pro_jump;        
        if isdebug; disp('Use Global_Var_RWRM'); end 
    end
    if ~exist('tau_vec','var') || isempty(tau_vec)
        tau_vec = repmat(1/L_net,L_net,1); 
    end        
    P0_L = sparse( kron( reshape(tau_vec,[],1), reshape(P0,[],1)));      
    Pt = A_RWRplus(Global_Var_RWRM.A_nLxnL_Normalized, restart, P0_L , [],[], IsNormalized);   
    if isa(AdjSet,'struct')
        AdjSetfieldnameset = fieldnames( AdjSet);  
    elseif isa(AdjSet,'table')
        AdjSetfieldnameset = AdjSet.Properties.VariableNames; 
    end 
    Pset =reshape(Pt,N_node, L_net) ; 
    TableScores = table; 
    if L_net>1     
        TableScores.(['ScoreGeoMean' ] ) = geomean( Pset , 2  ) ;     
    end
    for ii_L_net = 1: L_net
        TableScores.([AdjSetfieldnameset{ii_L_net}] ) = Pset(:, ii_L_net) ;    
    end   
    
end
