classdef (Abstract) SpinGroup < handle
    %SPINGROUP The parent class for all spin groups
    %   defines the interface of a spin group
    
    properties (Abstract)
    end
    
    methods (Abstract)
        dm = dM(self,position,M,time,PS,B0)
        % dM(self,position,M,time) - calulates the dm of the spin at some
        % position and time in the calculation frame defined by the spin.
        val = calculationFrame(self,B0)
        % calculationFrame the vector defining the angular momentum of the
        % calculation frame dm is defined in. remember to create a get method
        % for this variable in in inherited classes
    end
    
end

