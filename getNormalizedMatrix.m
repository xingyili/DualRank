% getNormalizedMatrix() -  get normalized matrix
%
% Usage:
%              >> getNormalizedMatrix(Adj, NormalizationType, SetIsolatedNodeSelfLoop );
%              >> WAdj = getNormalizedMatrix(Adj, NormalizationType, SetIsolatedNodeSelfLoop )
% Inputs:
%       Adj   = adjacency matrix of networks
%       NormalizationType   = normalization type
%       SetIsolatedNodeSelfLoop   = set isolated node
% Outputs:
%       WAdj  = normalized matrix
function WAdj = getNormalizedMatrix(Adj, NormalizationType, SetIsolatedNodeSelfLoop )
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %   
    if ischar(NormalizationType)
    %         NormalizationType =  (NormalizationType);
        switch NormalizationType 
            case { 'row' ,'ProbabilityNormalizationRow',    }
                NormalizationName = 'ProbabilityNormalization' ;  %  'Random Walk'  
                dim =2; 
            otherwise
                error(['There is no type of normalization: ',char( string(NormalizationType) )] );
        end
    elseif isnumeric(  NormalizationType   ) 
        NormalizationName = 'ProbabilityNormalization' ;  %  'Random Walk'  
        dim = NormalizationType; 
    end 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     %
    switch NormalizationName 
        case 'ProbabilityNormalization' 
            degrees = sum(Adj,dim);
            if any( degrees~=1)
                WAdj = Adj./ ( degrees+eps  );           
                % % WAdj = Adj./ repmat( degrees +eps,[size(Adj,1),1]); 
            else
                WAdj = Adj; 
            end
            % 
            if SetIsolatedNodeSelfLoop
                ii = find( ~degrees ); 
                idx = sub2ind( size(Adj), ii,ii ); 
                WAdj(idx) = 1;  % set to be 1 for isolated nodes, 
            end 
        otherwise
            error(['NormalizationName is wrong: ',char(string(NormalizationName) )   ]);
    end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
 
end
