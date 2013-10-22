% M_Map - mapping toolbox (Author: rich@ocgy.ubc.ca)
% Version 1.0 7/May/1997
%
% You have collected your data, loaded it into Matlab, analyzed 
% everything to death, and now you want to make a simple map showing 
% how it relates to the world. 
%
% But you can't. 
%
% Instead you have to figure out how to save all your data, and 
% then read it into a mapping program, and then spend all that extra 
% time figuring out why the mapping program doesn't give you what
% you expected it would...
%
% No more! 
%
%                            Announcing M_Map v1.0! 
%
% M_Map is a set of mapping tools written for Matlab v5. These include: 
%
%    1. Routines to project data in 13 different spherical 
%       projections (and determine inverse mappings) 
%    2. A grid generation routine to make nice axes with 
%       limits either in long/lat terms or in planar
%       X/Y terms. 
%    3. A coastline database (with 1/4 degree resolution) 
%    4. A global elevation database (1 degree resolution) 
%
%
% M_Map v1.0 is available via the web at 
%
%       http://www.ocgy.ubc.ca/~rich/map.html
%
%
% Toolbox contents
%
%    Contents.m - This file
%    m_demo.m   - demonstrates a few different maps.
%
%  User-callable functions
%
%    m_proj.m   - initializes projections
%    m_grid.m   - draws grids
% 
%    m_ungrid.m - erases grids (if you want to change grid parameters)
%
%    m_coast.m  - draws a coastline
%    m_elev.m   - draws elevation data
%
%    m_ll2xy.m  - converts from long/lat to map coordinates
%    m_xy2ll.m  - converts from map coordinates to long/lat
%
%  Internal functions (not meant to be user-callable)
%
%    private/mp_azim.m  - azimuthal projections
%    private/mp_cyl.m   - cylindrical projections (equatorial)
%    private/mp_conic.m - conic projections
%    private/mp_tmerc.m - transverse cylindrical projections
%    private/mp_omerc.m - oblique cylindrical projection
%
%    private/mu_util.m  - various utility routines
%
%    private/contourf.m - patched version of contourf 
%                         (matlab v5.0 version works incorrectly)    
%
%    private/m_coasts.mat - coastline data
%
%  HTML documentation
%
%    map.html           - Home page, examples
%    private/mapug.html - User's guide
%    private/*gif       - examples.
%  
%
% Questions or problems; email me - rich@ocgy.ubc.ca.
%
% Rich Pawlowicz
% Oceanography, Dept. of Earth and Ocean Sciences, Univ. of British Columbia, 
% 6270 University Blvd., Vancouver, B.C. CANADA V6T 1Z4
% email: rich@ocgy.ubc.ca 


    

