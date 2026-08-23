function x = elm_get_var(fname,varname,levelindex,timeindex,varargin)
%% DART elm_get_var - reconstitute an ELM vector field on its lon/lat grid.
%
% EXAMPLE 1:
% fname      = 'elm_restart.nc';
% varname    = 'T_SOISNO';
% levelindex = 1;
% timeindex  = 1;
% x          = elm_get_var(fname,varname,levelindex,timeindex);
% Optional name/value argument:
% x = elm_get_var(fname,varname,levelindex,timeindex,'verbose','no');
% % 'frac_sno',    'QTY_SNOWCOVER_FRAC',
% % 'DZSNO',       'QTY_SNOW_THICKNESS',
% % 'H2OSOI_LIQ',  'QTY_LIQUID_WATER',
% % 'H2OSOI_ICE',  'QTY_ICE',
% % 'T_SOISNO',    'QTY_SOIL_TEMPERATURE',
%
% 'frac_sno',    'QTY_SNOWCOVER_FRAC',
% 'DZSNO',       'QTY_SNOW_THICKNESS',
% 'H2OSOI_LIQ',  'QTY_LIQUID_WATER',
% 'H2OSOI_ICE',  'QTY_ICE',
% 'T_SOISNO',    'QTY_SOIL_TEMPERATURE',

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

p = inputParser;

addRequired(p,'fname',      @ischar);
addRequired(p,'varname',    @ischar);
addRequired(p,'levelindex', @isnumeric);
addRequired(p,'timeindex',  @isnumeric);

default_verbosity = 'yes';

if verLessThan('matlab','R2013b')
    addParamValue(p,'verbose', default_verbosity, @ischar);    %#ok<NVREPL>
   else
    addParameter(p,'verbose',  default_verbosity, @ischar);
   end
p.parse(fname, varname, levelindex, timeindex, varargin{:});

if ~isempty(fieldnames(p.Unmatched))
    disp('Extra inputs:')
    disp(p.Unmatched)
end

global verbose
if (strncmpi(p.Results.verbose,'y',1))
    verbose = 1;
else
    verbose = 0;
end

if (exist(p.Results.fname,'file') ~= 2)
    error('%s does not exist.',p.Results.fname)
end

if (~nc_isvar(p.Results.fname,p.Results.varname))
   error('%s does not have a "%s" variable.',fname,varname)
end

x.filename   = p.Results.fname;
x.varname    = p.Results.varname;
x.levelindex = p.Results.levelindex;
x.timeindex  = p.Results.timeindex;


%% define and get the basic coordinate variables

%variables = {'area', 'lon', 'lat', 'landfrac'};
%
%for ivar = 1:length(variables)
%   if (~nc_isvar(fname,variables{ivar}))
%      error('%s does not have a "%s" variable.',fname,variables{ivar})
%   else
%      expr = sprintf('x.%s = ncread(fname,variables{ivar});',variables{ivar});
%      eval(expr)
%   end
%end

%% Determine whether the variable is partitioned by gridcell, landunit, column, or PFT.

% The weights sum to 1.0, so no further normalization is needed.

vinfo = ncinfo(fname, x.varname);
ndims = length(vinfo.Dimensions);
weights = [];
for idim = 1:ndims
    switch lower(vinfo.Dimensions(idim).Name)
      case {'time','copy','levgrnd','levsoi','levdcmp', ...
            'levlak','levtot','levsno','levsno1'}
      case {'gridcell'}
            weights  = ones(size(ncread(fname,'grid1d_ixy')));
            loninds  = ncread(fname,'grid1d_ixy');
            latinds  = ncread(fname,'grid1d_jxy');
            gridlons = ncread(fname,'lon');
            gridlats = ncread(fname,'lat');
            lons     = gridlons(loninds);
            lats     = gridlats(latinds);
            ambigdim = idim; %#ok<NASGU>
      case {'landunit'}
            weights  = ncread(fname,'land1d_wtxy');
            loninds  = ncread(fname,'land1d_ixy');
            latinds  = ncread(fname,'land1d_jxy');
            gridlons = ncread(fname,'lon');
            gridlats = ncread(fname,'lat');
            lons     = gridlons(loninds);
            lats     = gridlats(latinds);
            ambigdim = idim; %#ok<NASGU>
      case {'column'}
            weights  = ncread(fname,'cols1d_wtxy');
            loninds  = ncread(fname,'cols1d_ixy');
            latinds  = ncread(fname,'cols1d_jxy');
            lons     = ncread(fname,'cols1d_lon');
            lats     = ncread(fname,'cols1d_lat');
 %          snlsno   = ncread(fname,'SNLSNO'); % fill value is -9999
            ambigdim = idim; %#ok<NASGU>
      case {'pft'}
            if nc_isvar(fname,'pfts1d_wtgcell')
                weights  = ncread(fname,'pfts1d_wtgcell');
            elseif nc_isvar(fname,'pfts1d_wtxy')
                weights  = ncread(fname,'pfts1d_wtxy');
            else
                error('unable to determine PFT weights for "%s"',x.varname)
            end
            loninds  = ncread(fname,'pfts1d_ixy');
            latinds  = ncread(fname,'pfts1d_jxy');
            lons     = ncread(fname,'pfts1d_lon');
            lats     = ncread(fname,'pfts1d_lat');
            ambigdim = idim; %#ok<NASGU>
        case {'lon','lat'}
            disp('Tim ... fix this.')
            error('unsupported gridded data "%s"',vinfo.Dimensions(idim).Name)
      otherwise
            error('unsupported dimension "%s"',vinfo.Dimensions(idim).Name)
   end
end

if ( exist('ambigdim','var') ~= 1)
   error('Could not determine weights/indexes for %s',varname)
end

%% unpack the sparse representation into a matrix.

myinfo.fname      = fname;
myinfo.timeindex  = timeindex;
myinfo.levelindex = levelindex;
[start,count]     = GetNCindices(myinfo,'fname',varname);
x.dat             = ncread(fname,varname,start,count);

% If the variable is a snow variable, use the information in snlsno
% when they erode values, it is not clear that they reset them to the
% _FillValue

%% unpack the sparse representation into a matrix.
x.datmat   = zeros(max(latinds),max(loninds));
x.wgtmat   = zeros(size(x.datmat));
x.lonarray = zeros(max(loninds),1);
x.latarray = zeros(max(latinds),1);

for i = 1:length(x.dat)
   ilon = loninds(i);
   jlat = latinds(i);

   if (isfinite(x.dat(i)))
      x.datmat(jlat,ilon) = x.datmat(jlat,ilon) + weights(i)*x.dat(i);
      x.wgtmat(jlat,ilon) = x.wgtmat(jlat,ilon) + weights(i);
        x.lonarray(ilon) = lons(i);
        x.latarray(jlat) = lats(i);
   end
end

empty_cells           = find(x.wgtmat == 0.0);
x.datmat(empty_cells) = NaN;

% =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

function yesno = nc_isvar(fname,varname)

finfo = ncinfo(fname);

nvars = length(finfo.Variables);

yesno = -1;  % variable is not present

for ivar = 1:nvars

    yesno = strcmp(varname, finfo.Variables(ivar).Name);

    if yesno
        return
   end

end
