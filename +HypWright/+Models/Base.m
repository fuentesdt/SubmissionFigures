classdef (Abstract) Base
    %BASE the interface class for all models
    %   This is the abstract class that all models will inheret from.
    
    properties (Abstract)
    end
    methods (Abstract,Static)
        evaluate()
        fitData()
    end
    methods (Abstract,Static, Access = protected)
        parseParams()
        fitFunction()
    end
end

