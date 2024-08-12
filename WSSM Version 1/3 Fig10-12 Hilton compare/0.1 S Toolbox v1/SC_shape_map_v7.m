function g = SC_shape_map_v7(input)
% Loads the map z=g(zeta) from the database -- see corresponding folder.
% SHAPE DIRECTORY
%   11: square
%   12: L-shape
%   13: three-lemniscate
%
%   The following correspond to figure numbers from [1], e.g. 12 = figure 1b.
%   92: rounded snowflake
%   93: smooth random boundary
%   101: bean
%   102: blade
%
%   The following correspond to wildfire data from previous works.
%   12144: Viegas et al 2012 Fig 14d
%   173: Sharples and Hilton 2017 Fig 3
%   185: Hilton et al 2018 Fig 5
%   187: Hilton et al 2018 Fig 5
%
% REFERENCES
%   [1] Gopal, A., & Trefethen, L. N. (2019). Representation of conformal 
%       maps by rational functions. Numerische Mathematik, 142, 359-382.
%   (See the references therein for the relevant shapes).
%
% END OF DOCUMENTATION
%
%Code
% --- Wildfire shapes ---
if input==185 %Hilton18-fig5 fire
    load('g_hilton18_fig5b_v7.mat','g');

elseif input==187 %Hilton18-fig7 fire
    load('g_hilton18_fig7e_csi_v7.mat','g'); 
end
end