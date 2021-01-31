% getMultiplexMatrixFromAdjSet() -  Get multiplex matrix from the adjacency
%
% Usage:
%              >> getMultiplexMatrixFromAdjSet(AdjSet, pro_jump, NormalizationType, InterLayerMappingType);
%              >> [A_nLxnL,N_node, L_net] = ...
%                       getMultiplexMatrixFromAdjSet(AdjSet, pro_jump, NormalizationType, InterLayerMappingType);
% Inputs:
%       AdjSet   = adjacency matrix of networks
%       pro_jump   =  the parameter that denotes the relative importance of the states of the neighbors in the layer or the states of nodes in other layers
%       NormalizationType   = normalization type
%       InterLayerMappingType   = inter layer mapping type
% Outputs:
%       A_nLxnL, A_nLxnL, L_net  = multiplex network information
function [A_nLxnL,N_node, L_net] = getMultiplexMatrixFromAdjSet(AdjSet, pro_jump, NormalizationType, InterLayerMappingType)
    if ~exist('NormalizationType','var') 
        NormalizationType = 'col';
    elseif isempty(NormalizationType) || strcmpi(NormalizationType,'none') 
        NormalizationType = 'None';
    end
    %
    if ~exist('InterLayerMappingType','var') 
        InterLayerMappingType ='FullMapping';
    else
        InterLayerMappingType ='PartialMapping'; % 
    end
%
    delta =pro_jump; 
    SetIsolatedNodeSelfLoop = true; 
    if isa(AdjSet,'struct') || isa(AdjSet,'table')
        if isa(AdjSet,'struct')
            fieldnameset = fieldnames( AdjSet); 
        elseif isa(AdjSet,'table')
            fieldnameset = AdjSet.Properties.VariableNames; 
        end
        L_net  =  length( fieldnameset ); 
        N_node =  length( AdjSet.(fieldnameset{1}));
        NxL = N_node*L_net ; 
        %
        if L_net==1;  delta = 0; end %% no jumping for only single layer. 
        %        
        if strcmpi(InterLayerMappingType, 'PartialMapping' )
            k_colvec = []; k_rowvec = [];                
            for ii_net = 1: L_net
                    k_colvec = [k_colvec; sum(AdjSet.(fieldnameset{ii_net}),2) ] ;    
                    k_rowvec = [k_rowvec, sum(AdjSet.(fieldnameset{ii_net}),1) ] ;          
            end  
            A_nLxnL = repmat( speye( N_node, N_node), L_net,L_net); 
            %
            Kcutinter = 0 ; 
            A_nLxnL(k_colvec<=Kcutinter,:  ) = 0;  
            A_nLxnL(: , k_rowvec<=Kcutinter) = 0; 
            %
            ind = sub2ind(  size(A_nLxnL) , [1:NxL],[1:NxL] );  
            A_nLxnL(ind) =0;   %
            switch NormalizationType 
                case { 'row' ,'ProbabilityNormalizationRow',    }
                    k_colvec_interlayer = sum( A_nLxnL , 2 );   
                    A_nLxnL = A_nLxnL.*( delta./sparse(k_colvec_interlayer +eps) ) ;   
                    ind_setselfloop = k_colvec==0;  
                    
                case 'None'   % no operation  
                    A_nLxnL = delta.*A_nLxnL; 
                    ind_setselfloop = k_colvec==0 & (k_rowvec'==0);  

                otherwise
                    error(['There is no type of normalization: ',char( string(NormalizationType) )] );
            end   
            % 
            for ii_net = 1: L_net
                idx = N_node*(ii_net-1)+[1: N_node ] ; 
                if strcmpi(NormalizationType,'None') 
                    A_nLxnL(idx,idx)=  (1-delta).*AdjSet.(fieldnameset{ii_net}); 
                else
                    A_nLxnL(idx,idx)=  (1-delta).*getNormalizedMatrix( AdjSet.(fieldnameset{ii_net}) , NormalizationType, false );  % 对孤立节点 下面单独设置 
                end                   
            end
            if SetIsolatedNodeSelfLoop
                ii = find( ind_setselfloop ); 
                ind = sub2ind( size(A_nLxnL), ii,ii ); 
                A_nLxnL(ind) = 1;  % set to be 1 for isolated nodes, 
            end        
            
        elseif strcmpi(InterLayerMappingType, 'FullMapping' )
            A_nLxnL = repmat( (delta/(L_net-1+eps)).*speye( N_node, N_node), L_net,L_net);            
            for ii_net = 1: L_net
                idx = N_node*(ii_net-1)+[1: N_node ] ; 
                if strcmpi(NormalizationType,'None') 
                    A_nLxnL(idx,idx)=  (1-delta).*AdjSet.(fieldnameset{ii_net}); 
                else
                    A_nLxnL(idx,idx)=  (1-delta).*getNormalizedMatrix( AdjSet.(fieldnameset{ii_net}) , NormalizationType, SetIsolatedNodeSelfLoop ); 
                end                   
            end 
        else
            error(['There is no type of InterLayerType: ',char( string(InterLayerMappingType) )] );
        end        
    elseif iscell( AdjSet )
        L_net  = numel( AdjSet); 
        N_node = length( AdjSet{1});         
        if L_net==1;  delta = 0; end %% no jumping for only single layer. 

        switch NormalizationType 
            case { 'row' ,'ProbabilityNormalizationRow',    }
                A_nLxnL = repmat( (delta/(L_net-1+eps)).*speye( N_node, N_node), L_net,L_net);   
                
            case 'None'
                A_nLxnL = repmat( ( delta ).*speye( N_node, N_node), L_net,L_net);  
                
            otherwise
                error(['There is no type of normalization: ',char( string(NormalizationType) )] );
        end  
        
        for ii_net = 1: L_net
            idx = N_node*(ii_net-1)+[1: N_node ] ; 
            if strcmpi(NormalizationType,'None') 
                A_nLxnL(idx,idx)=  (1-delta).*AdjSet{ii_net} ;   
            else
                A_nLxnL(idx,idx)=  (1-delta).*getNormalizedMatrix( AdjSet{ii_net} , NormalizationType, SetIsolatedNodeSelfLoop );  
            end 
        end 
        
    else
        error(['AdjSet is wrong. It should be a cell matrix or struct.' ]);
        
    end 
        
end
