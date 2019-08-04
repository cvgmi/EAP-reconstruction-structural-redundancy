classdef surflet
    properties
        adjoint = 0;
        Vx=0; Vy = 0; Vz = 0;
        bo = 0;
        Pyr_mode = 0;
        Lev_array={};
        Y = 0;
        Recinfo = 0;
        size_info=0;
        
    end
    
    methods
        function obj = surflet(volume) %constructor
            Level_64 = [-1 3 3; 3 -1 3; 3 3 -1]; % 3 * 64 directions
            Level_16 = [-1 2 2; 2 -1 2; 2 2 -1]; % 3 * 16 directions
            Level_4 =  [-1 1 1; 1 -1 1; 1 1 -1]; % 3 * 4 directions
            Level_1 =  [-1 0 0; 0 -1 0; 0 0 -1]; % 3 directions (i.e. the hourglass decomposition)
            
            [obj.Vx, obj.Vy , obj.Vz] = size(volume);
            obj.bo = 12;
            if length(volume)==16
                obj.Pyr_mode = 1;
                obj.Lev_array = {Level_1,Level_1};
            elseif length(volume)==32
                obj.Pyr_mode = 1;
                obj.Lev_array = {Level_1,Level_1};
            elseif length(volume)==64
                obj.Pyr_mode = 2;
                obj.Lev_array = {Level_1,Level_1};
            elseif length(volume)==128
                obj.Pyr_mode = 2;
                obj.Lev_array = {Level_4, Level_16};
            elseif length(volume)==256
                obj.Pyr_mode = 1.5;
                obj.Lev_array = {Level_4, Level_16};
            else
                display('Unknown size...');
                return;
            end
            
            X0 = zeros(size(volume));
            
            [tmp, obj.Recinfo] = Surfdec(volume, obj.Pyr_mode, obj.Lev_array, 'ritf', 'bo', obj.bo);
            [~, obj.size_info] = surf_coeff2vec(tmp);
            obj.adjoint = 0;
        end
       
    end
end


