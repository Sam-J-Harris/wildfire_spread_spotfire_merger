function g = shape_SC_v1(input)
% Loads the map z=g(zeta) from the database -- see corresponding folder.
% SHAPE DIRECTORY
%   11: square
%   12: L-shape
%
%   The following correspond to figure numbers from [1], e.g. 12 = figure 1b.
%   92: rounded snowflake
%   101: bean
%   102: blade
%
%   The following correspond to wildfire data from [2].
%   185: [2] Hilton et al 2018 Fig 5
%   187: [2] Hilton et al 2018 Fig 7
%
% REFERENCES
%   [1] Gopal, A., & Trefethen, L. N. (2019). Representation of conformal 
%       maps by rational functions. Numerische Mathematik, 142, 359-382.
%   (See the references therein for the relevant shapes).
%
% END OF DOCUMENTATION
%
%Code
if input==11 %square
    load('g_square.mat','g'); 

elseif input==12 %L-shape
    load('g_lshape.mat','g'); 

% --- Shapes from Gopal and Trefethen 2019 ---
elseif input == 92 %rounded snowflake
    load('g_rounded_snowflake.mat','g');  

elseif input==101 %bean
    load('g_bean.mat','g'); 

elseif input==102 %blade
    load('g_blade.mat','g');  

% --- Wildfire shapes ---
elseif input==185 %Hilton18-fig5 fire
    load('g_hilton18_fig5.mat','g');

elseif input==187 %Hilton18-fig7 fire
    load('g_hilton18_fig7.mat','g'); 
end
end