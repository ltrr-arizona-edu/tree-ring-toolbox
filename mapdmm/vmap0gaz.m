function varargout = vmap0gaz(varargin)

%VMAP0GAZ searches for named entries in the Vector Map Level 0 
%
%   gstruct = VMAP0GAZ('library','theme','topolevel','object') searches  
%   the Vector Map Level 0 CD-ROM for items which have names or strings
%   beginning with the object string.  There are four CD's, one for each
%   library: 'NOAMER' (North America), 'SASAUS' (Southern Asia and
%   Australia), 'EURNASIA' (Europe and Northern Asia) and 'SOAMAFR' (South
%   America and Africa).  The desired theme is specified by a two letter 
%   code string.  A list of valid codes will be displayed when an invalid 
%   code, such as '?', is entered. Topolevel may be one of the following:
%   'path', 'line', 'point' or 'text'. Matching items are returned in a
%   geographic data structure. Among the fields of this structure are the
%   latitude and longitude vectors associated with the items, and the
%   scalar mean latitude and longitude.
%  
%   gstruct = VMAP0GAZ('devicename','library','theme','topolevel','object') 
%   specifies the logical device name of the CD-ROM for computers which 
%  require them.
%   
% See also VMAP0DATA, VMAP0READ, DCWDATA, DCWGAZ, MLAYERS, DISPLAYM, EXTRACTM

%  Copyright (c) 1996-98 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  W. Stumpf
%   $Revision: 1.6 $    $Date: 1997/03/28 23:11:49 $

% check for valid inputs

% devicename, library, theme,object

if nargin == 4
	library	= varargin{1};
	theme 	= varargin{2};
	topolevel= varargin{3};
	object	= varargin{4};
	devicename = 'VMAP';
elseif nargin == 5
	devicename = varargin{1};
	library = varargin{2};
	theme = varargin{3};
	topolevel = varargin{4};
	object = varargin{4};
else
	error('Incorrect number of input arguments')
end


%  Check for valid topology level inputs

if isstr(topolevel) & strcmp(topolevel,'all'); topolevel = {'all'}; end

if isstr(topolevel)

	topoindx = strmatch(topolevel,{'patch','line','point','text'},'exact');
	if isempty(topoindx); error('Valid topology levels are ''patch'', ''line'', ''point'' or ''text'' ');end
	if nargout ~= 1
		error('Number of output arguments doesn''t match requested topology levels')
	end

elseif iscell(topolevel)

	if length(topolevel)==1 & strcmp(topolevel,'all');
		topolevel = {'patch','line','point','text'};
	end

	topoindx = [];

	for i=1:length(topolevel)
		if ~isstr(topolevel{i}); error('Topology level must be a cellarray of strings'); end
		thistopoindx = strmatch(topolevel{i},{'patch','line','point','text'},'exact');
		if isempty(thistopoindx); error('Topology level must be ''patch'', ''line'', ''point'' or ''text'' ');end
		topoindx = [topoindx thistopoindx];
	end % for

	topoindx = unique(topoindx);
	if length(topoindx) ~=length(topolevel)
		error('Redundent requests in topology levels')
	elseif nargout ~= length(topoindx)
		error('Number of outputs doesn''t match requested topology levels')
	end
else
	error('Topology level must be a string or cellarray of strings')
end % if


% build the pathname so that [pathname filename] is the full filename

filepath = fullfile(devicename,'VMAPLV0',library,filesep);
if isunix; 
	filepath = fullfile(devicename,'vmaplv0',lower(library),filesep); 
end

% check that we can get into disk

dirstruc = dir(filepath);
if length(dirstruc) == 0; error(['Disk ' library ' not mounted']); end

% check for valid theme request

CAT = vmap0read(filepath,'CAT');

if isempty(strmatch(theme,{CAT.coverage_name},'exact'))

	linebreak = double(sprintf('\n'));
	goodthemes = [strvcat(CAT.coverage_name) ...
				char(58*ones(length(CAT),1)) ...
			  	char(32*ones(length(CAT),1)) ...
				strvcat(CAT.description) ...
				char(linebreak*ones(length(CAT),1)) ];
	goodthemes = goodthemes';
	goodthemes = goodthemes(:);
	goodthemes = goodthemes';

	error(['Theme not present in library ' library char(linebreak) char(linebreak) ...
			 'Valid theme identifiers are: ' char(linebreak) ...
			 goodthemes ...
			])
end

if strcmp(lower(library),'rference')

% BROWSE layer is untiled

	dotiles = 1;
	tFT(1).tile_name = '';
	tFT(1).fac_id = 1;
else

% Get the essential libref information (tile name/number, bounding boxes)

	filepath = fullfile(devicename,'VMAPLV0',library,'TILEREF',filesep);
	if isunix; 
		filepath = fullfile(devicename,'vmaplv0',lower(library),'tileref',filesep); 
	end

	tFT = vmap0read(filepath,'TILEREF.AFT');

end

% Here is where the value description and feature tables reside

themepath = fullfile(devicename,'VMAPLV0',library,theme,filesep);
if isunix; 
	themepath = fullfile(devicename,'vmaplv0',lower(library),theme,filesep); 
end

% get a list of files in the requested directory

dirstruc = dir(themepath);
names = {dirstruc.name};
names = uppercell(names); % because pc converts to lowercase

% loop over points, lines, text
% handle faces separately

savestruct = [];

VDT = [];
cVDT = [];

% loop over feature topology level

for j=1:4

% construct the feature table name. We know it will be of form
% '**'POINT.PFT', etc.

	switch j

% Patch
	
	case 1
		if any(topoindx ==1)
			[savestruct,VDT,cVDT] = facegaz(object,savestruct,themepath,tFT,VDT,cVDT);	
		end
				
% Line		
	case 2

		if any(topoindx ==2)
			[savestruct,VDT,cVDT] = edgegaz(object,savestruct,themepath,tFT,VDT,cVDT);	
		end
		
% Point

	case 3


		if any(topoindx ==3)

			[savestruct,VDT,cVDT] = pointgaz(object,savestruct,themepath,tFT,VDT,cVDT);	
				
		end % if
% Text			
	otherwise

		if any(topoindx ==4)

			[savestruct,VDT,cVDT] = textgaz(object,savestruct,themepath,tFT,VDT,cVDT);						
				
		end % if topoindx

	end % switch

end % for j

varargout{1} = savestruct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function carray = uppercell(carray)
%UPPERCELL converts a cellarray of strings to uppercase

for i=1:length(carray);carray{i} = upper(carray{i});end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function carray = lowercell(carray)
%LOWERCELL converts a cellarray of strings to lowercase

for i=1:length(carray);carray{i} = lower(carray{i});end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [savelat,savelong] = vmap0gbfac(facenum,EDG,FAC,RNG)
%VMAP0GFAC builds a Digital Chart of the World Browse layer face from its edges

if facenum ==1;

% outer ring of the universe face cannot be displayed

	savelat = [];
	savelong = [];
	return
end

ringptrs = find([RNG.face_id]==facenum);

savelat = [];
savelong = [];
for i=1:length(ringptrs)

	checkface = RNG(ringptrs(i)).face_id;
	if checkface ~= facenum; warning('Face and Ring tables inconsistent');return; end

	startedge = RNG(ringptrs(i)).start_edge;

	[lat,long] = vmap0gbrng(EDG,startedge,facenum);


	if isempty(savelat)
		savelat = lat;
		savelong = long;
	else
		savelat = [savelat;NaN;lat];
		savelong = [savelong;NaN;long];
	end

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [savelat,savelong] = vmap0gbrng(EDG,startedge,facenum)

if nargin ~= 3
        error('Incorrect number of arguments')
end


% Follow face until we end up where we started. If left and right face
% are identical, use left/right edges to break the tie; if that's not
% unambiguous, use start/end nodes

llchunk = EDG(startedge).coordinates;

if  EDG(startedge).right_face(1) == facenum
        nextedge = EDG(startedge).right_edge(1);  % stay on this tile
        savelat = llchunk(:,2);
        savelong = llchunk(:,1);
elseif EDG(startedge).left_face(1) == facenum
        nextedge = EDG(startedge).left_edge(1);
        savelat  = flipud(llchunk(:,2));
        savelong = flipud(llchunk(:,1));
else
        warning('Error following face.');return
end % if

lastedge=startedge;
lastnodes = [EDG(startedge).start_node(1)  ...
                                EDG(startedge).end_node(1)  ];

while ~(nextedge == startedge)

        curnodes = [EDG(nextedge).start_node EDG(nextedge).end_node];

        llchunk = EDG(nextedge).coordinates;
        lat=llchunk(:,2);
        long=llchunk(:,1);

        rface = EDG(nextedge).right_face(1);
        lface = EDG(nextedge).left_face(1);

        if lface == rface

% Breaking tie with nextedge

                if EDG(nextedge).right_edge(1) == nextedge

                        savelat = [savelat;lat;flipud(lat)];
                        savelong = [savelong;long;flipud(long)];
                        lastedge = nextedge;
                        nextedge = EDG(nextedge).left_edge(1);

                elseif EDG(nextedge).left_edge(1) == nextedge

                        savelat = [savelat;flipud(lat);lat];
                        savelong = [savelong;flipud(long);long];
                        lastedge = nextedge;
                        nextedge = EDG(nextedge).right_edge(1);  %stay on this tile

                else

% Breaking tie with nodes

                        starttest = find(curnodes(1) == lastnodes);
                        endtest = find(curnodes(2) == lastnodes);

                        if isempty(starttest) & ~isempty(endtest)

                                savelat = [savelat;flipud(lat)];
                                savelong = [savelong;flipud(long)];
                                nextedge = EDG(nextedge).left_edge(1);

                        elseif ~isempty(starttest) & isempty(endtest)

                                savelat = [savelat;lat];
                                savelong = [savelong;long];
                                nextedge = EDG(nextedge).right_edge(1);  % stay on this tile

                        else
                                warning('Error following face..')
								return

                        end % if
                end % if

        elseif rface == facenum
                savelat = [savelat;lat];
                savelong = [savelong;long];
                lastedge = nextedge;
                nextedge = EDG(nextedge).right_edge(1);  % stay on this tile
        elseif lface == facenum
                savelat = [savelat;flipud(lat)];
                savelong = [savelong;flipud(long)];
                lastedge = nextedge;
                nextedge = EDG(nextedge).left_edge(1);
        else
                warning('Error following face...')
                return
        end % if

        nextedge = nextedge(1);
        lastnodes = curnodes;

end % while

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [description,VDT,cVDT] = descript(entry,FT,FTfield,FTfilename,themepath,VDT,cVDT)

% build up description codes

description = '';

for i=1:length(FTfield)

	fieldname = FTfield(i).name;
	
	switch FTfield(i).VDTname
	
	case 'char.vdt'
		
		if isempty(cVDT)
			dirstruc = dir(themepath);
			names = {dirstruc.name};
			names = uppercell(names); % because pc converts to lowercase
			if any(strcmp(names,'CHAR.VDT'))==1
				cVDT = vmap0read(themepath,'CHAR.VDT');
			end
		end
	
		fieldname = FTfield(i).name;

		strng = getfield(FT,{entry},fieldname);
		
		if ~isempty(cVDT)
			tableindx = strmatch(lower(FTfilename), {cVDT.table});
			fieldindx = strmatch(fieldname, {cVDT.attribute});
			valueindx = strmatch(strng,{cVDT.value});
			cVDTindx = intersect(tableindx, intersect(fieldindx,valueindx));
			
			if ~isempty(cVDTindx)
				strng = cVDT(cVDTindx(1)).description;
			end
		end
		
		description = [description '; ' deblank(strng)];


	case 'int.vdt'

		if isempty(VDT)
			dirstruc = dir(themepath);
			names = {dirstruc.name};
			names = uppercell(names); % because pc converts to lowercase
			if any(strcmp(names,'INT.VDT'))==1
				VDT = vmap0read(themepath,'INT.VDT');
			end
		end
			
		val = getfield(FT,{entry},fieldname);
		strng = num2str(val);

		tableindx = strmatch(lower(FTfilename), {VDT.table});
		fieldindx = strmatch(fieldname, {VDT.attribute});
		valueindx = find(val == [VDT.value]);
		VDTindx = intersect(tableindx, intersect(fieldindx,valueindx));
		
		if ~isempty(VDTindx)
			strng = VDT(VDTindx(1)).description;
			if strmatch(strng,'Unknown')  % expand labels of things with value unknown, 
				strng = [ FTfield(i).description ': ' strng];
			end
		end
	
		description = [description '; ' deblank(strng)];
	end % switch

end % for
		

% % extract value fields
% 
% 	if ~isempty(valfieldIDs) % have contour values 
% 
% 		for k = 1:length(valfieldIDs)
% 
% 			eval(['val = FT(indx(1)).' ftfields{valfieldIDs(k)}	';']);
% 			description = [description '; ' num2str(val)];
% 
% 		end %k
% 
% 	end % ~isempty(valfieldIDs)
% 
% 


description = description(3:length(description));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [savestruct,VDT,cVDT] = textgaz(object,savestruct,themepath,tFT,VDT,cVDT)


wildcard = '*.TFT';
if strcmp(computer,'PCWIN')| isunix; wildcard = lower(wildcard); end

FTfilenames = dir([themepath wildcard]);

% Loop on Area Feature Tables

if ~isempty(FTfilenames)

	for k=1:length(FTfilenames)

		FTfilename = FTfilenames(k).name;

% Read the feature table. The strings are at the entity tile level, so find matches there.
		
		[hitstruct,FTfield] = vmap0read(themepath,FTfilename);
		
% Having gathered all entries, go to the tiles and get the geographic data.
		
		if ~isempty(hitstruct)
		
			[savestruct,VDT,cVDT] = extracttext(object,savestruct,hitstruct,themepath,FTfield,FTfilename,tFT,VDT,cVDT);
		
		end % if ~isempty(hitstruct)
	
	end % do k				

end % if ~isempty(FTfilenames)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [savestruct,VDT,cVDT] = pointgaz(object,savestruct,themepath,tFT,VDT,cVDT)


wildcard = '*.PFT';
if strcmp(computer,'PCWIN')| isunix; wildcard = lower(wildcard); end

FTfilenames = dir([themepath wildcard]);

% Loop on Area Feature Tables

if ~isempty(FTfilenames)

	for k=1:length(FTfilenames)

		FTfilename = FTfilenames(k).name;

% warn of number of entries (for the moment)		
		
		xfilename = FTfilename;
		xfilename(length(xfilename))='X';
		if isunix; xfilename = lower(xfilename); end
		xfid = fopen([themepath xfilename],'r');
		if xfid ~= -1; 
			fclose(xfid);
			[X,nX] = vmap0rdx(themepath,xfilename);
			
			disp(['Reading feature table with ' num2str(nX) ' entries'])
		end
% Find enties of the feature table that match the desired object
		
		[hitstruct,FTfield] = findhits(object,themepath,FTfilename);
		
% Having gathered matching entries, go to the tiles and get the geographic data.
% Construct the patch, compute the mean position and append result to geographic
% data structure.
		
		if ~isempty(hitstruct)
		
			[savestruct,VDT,cVDT] = extractpoints(savestruct,hitstruct,themepath,FTfield,FTfilename,tFT,VDT,cVDT);
		
		end % if ~isempty(hitstruct)
	
	end % do k				

end % if ~isempty(FTfilenames)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [savestruct,VDT,cVDT] = edgegaz(object,savestruct,themepath,tFT,VDT,cVDT)


wildcard = '*.LFT';
if strcmp(computer,'PCWIN')| isunix; wildcard = lower(wildcard); end

FTfilenames = dir([themepath wildcard]);

% Loop on Area Feature Tables

if ~isempty(FTfilenames)

	for k=1:length(FTfilenames)

		FTfilename = FTfilenames(k).name;

% warn of number of entries (for the moment)		
		
		xfilename = FTfilename;
		xfilename(length(xfilename))='X';
		if isunix; xfilename = lower(xfilename); end
		xfid = fopen([themepath xfilename],'r');
		if xfid ~= -1; 
			fclose(xfid);
			[X,nX] = vmap0rdx(themepath,xfilename);
			
			disp(['Reading feature table with ' num2str(nX) ' entries'])
		end
	
% Find enties of the feature table that match the desired object
		
		[hitstruct,FTfield] = findhits(object,themepath,FTfilename);
		
% Having gathered matching entries, go to the tiles and get the geographic data.
% Construct the patch, compute the mean position and append result to geographic
% data structure.
		
		if ~isempty(hitstruct)
		
			[savestruct,VDT,cVDT] = extractedges(savestruct,hitstruct,themepath,FTfield,FTfilename,tFT,VDT,cVDT);
		
		end % if ~isempty(hitstruct)
	
	end % do k				

end % if ~isempty(FTfilenames)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [savestruct,VDT,cVDT] = facegaz(object,savestruct,themepath,tFT,VDT,cVDT)


wildcard = '*.AFT';
if strcmp(computer,'PCWIN')| isunix; wildcard = lower(wildcard); end

FTfilenames = dir([themepath wildcard]);

% Loop on Area Feature Tables

if ~isempty(FTfilenames)

	for k=1:length(FTfilenames)

		FTfilename = FTfilenames(k).name;

% warn of number of entries (for the moment)		
		
		xfilename = FTfilename;
		xfilename(length(xfilename))='X';
		if isunix; xfilename = lower(xfilename); end
		xfid = fopen([themepath xfilename],'r');
		if xfid ~= -1; 
			fclose(xfid);
			[X,nX] = vmap0rdx(themepath,xfilename);
			
			disp(['Reading feature table with ' num2str(nX) ' entries'])
		end
	
% Find enties of the feature table that match the desired object
		
		[hitstruct,FTfield] = findhits(object,themepath,FTfilename);
		
% Having gathered matching entries, go to the tiles and get the geographic data.
% Construct the patch, compute the mean position and append result to geographic
% data structure.
		
		if ~isempty(hitstruct)
		
			[savestruct,VDT,cVDT] = extractfaces(savestruct,hitstruct,themepath,FTfield,FTfilename,tFT,VDT,cVDT);
		
		end % if ~isempty(hitstruct)
	
	end % do k				

end % if ~isempty(FTfilenames)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hitstruct,FTfield] = findhits(object,themepath,FTfilename)

hitstruct = [];

% read first record from the feature table to see if it has a name field

[FT,FTfield] = vmap0read(themepath,FTfilename,1);
namefieldnum = strmatch('nam',{FTfield.name},'exact');
 
% if it does, read the full table, and find entries for which the 'nam' field 
% matches the desired name. Save those entries in a new table

if ~isempty(namefieldnum) % feature table has a name field
	FT = vmap0read(themepath,FTfilename);
	matchIDs = strmatch(lower(object),lowercell({FT.nam}));
	
% 			disp(strvcat(FT(matchIDs).nam))

	nhits = length(hitstruct);
	if isempty(hitstruct) & length(matchIDs) > 0; 
		hitstruct = FT(matchIDs(1)); % initialize to a structure so assignments possible
	end
	for i=1:length(matchIDs)
		hitstruct(nhits+i) = FT(matchIDs(i));
	end % for j
end % if
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [savestruct,VDT,cVDT] = extracttext(object,savestruct,hitstruct,themepath,FTfield,FTfilename,tFT,VDT,cVDT)

m= length(savestruct);
dotiles = unique([hitstruct.tile_id]');

for i=1:length(dotiles)
	
% Construct the path to the tile entity data

% 	tileid = find(dotiles(i)==[tFT.fac_id]);
	tileid = dotiles(i); % similar change to others?	

	tilepath = tFT(tileid).tile_name;
	tilepath = [ strrep(tilepath,'\',filesep) filesep];
	if isunix; tilepath = lower(tilepath); end
	tilepath = [themepath tilepath];
	tilepath = strrep(tilepath,[filesep filesep],filesep); % for the Browse library

% Read the entity table for that tile
				
	TXT = vmap0read(tilepath,'TXT');
	
% find matches

	hits = strmatch(object,lowercell({TXT.string}));

% 	hits = find([hitstruct.tile_id] == dotiles(i));

% Construct the patch for the hits on that tile
	
	for n=1:length(hits)
		indx = hits(n)
		lat = TXT(indx).shape_line(1,2);
		lon = TXT(indx).shape_line(1,1);
		
		[description,VDT,cVDT] = descript(hits(n),hitstruct,FTfield,FTfilename,themepath,VDT,cVDT);
		
		m=m+1;
		savestruct(m).string = TXT(indx).string;
		savestruct(m).lat = lat;
		savestruct(m).long = lon;
		savestruct(m).type = 'text';
		savestruct(m).tag = description;
		savestruct(m).altitude = [];
		savestruct(m).otherproperty = {};
	end

end % for i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [savestruct,VDT,cVDT] = extractpoints(savestruct,hitstruct,themepath,FTfield,FTfilename,tFT,VDT,cVDT)

m= length(savestruct);
dotiles = unique([hitstruct.tile_id]');

for i=1:length(dotiles)
	
% Construct the path to the tile entity data

	tileid = dotiles(i); % similar change to others?	
	hits = find([hitstruct.tile_id] == dotiles(i));
	
	tilepath = tFT(tileid).tile_name;
	tilepath = [ strrep(tilepath,'\',filesep) filesep];
	if isunix; tilepath = lower(tilepath); end
	tilepath = [themepath tilepath];
	tilepath = strrep(tilepath,[filesep filesep],filesep); % for the Browse library

% Read the face table for that tile
				
	END = vmap0read(tilepath,'END');

% Construct the patch for the hits on that tile
	
	for n=1:length(hits)
		indx = hitstruct(hits(n)).end_id ; % right?
		lat = END(indx).coordinate(:,2);
		lon = END(indx).coordinate(:,1);
		
		[latm,lonm] = meanm(lat(~isnan(lat)),lon(~isnan(lon)));
		
		[description,VDT,cVDT] = descript(hits(n),hitstruct,FTfield,FTfilename,themepath,VDT,cVDT);
		
		m=m+1;
		savestruct(m).string = hitstruct(hits(n)).nam;
		savestruct(m).lat = lat;
		savestruct(m).long = lon;
		savestruct(m).meanlat = latm;
		savestruct(m).meanlong = lonm;
		savestruct(m).type = 'line';
		savestruct(m).tag = description;
		savestruct(m).altitude = [];
		savestruct(m).otherproperty = {};
	end

end % for i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [savestruct,VDT,cVDT] = extractedges(savestruct,hitstruct,themepath,FTfield,FTfilename,tFT,VDT,cVDT)

m= length(savestruct);
dotiles = unique([hitstruct.tile_id]');

for i=1:length(dotiles)
	
% Construct the path to the tile entity data

	tileid = dotiles(i); % similar change to others?	
	hits = find([hitstruct.tile_id] == dotiles(i));
	
	tilepath = tFT(tileid).tile_name;
	tilepath = [ strrep(tilepath,'\',filesep) filesep];
	if isunix; tilepath = lower(tilepath); end
	tilepath = [themepath tilepath];
	tilepath = strrep(tilepath,[filesep filesep],filesep); % for the Browse library

% Read the face table for that tile
				
	EDG = vmap0read(tilepath,'EDG');

% Construct the patch for the hits on that tile
	
	for n=1:length(hits)
		indx = hitstruct(hits(n)).edg_id
		lat = EDG(indx).coordinates(:,2);
		lon = EDG(indx).coordinates(:,1);
		
		[latm,lonm] = meanm(lat(~isnan(lat)),lon(~isnan(lon)));
		
		[description,VDT,cVDT] = descript(hits(n),hitstruct,FTfield,FTfilename,themepath,VDT,cVDT);
		
		m=m+1;
		savestruct(m).string = hitstruct(hits(n)).nam;
		savestruct(m).lat = lat;
		savestruct(m).long = lon;
		savestruct(m).meanlat = latm;
		savestruct(m).meanlong = lonm;
		savestruct(m).type = 'line';
		savestruct(m).tag = description;
		savestruct(m).altitude = [];
		savestruct(m).otherproperty = {};
	end

end % for i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [savestruct,VDT,cVDT] = extractfaces(savestruct,hitstruct,themepath,FTfield,FTfilename,tFT,VDT,cVDT)

m= length(savestruct);
dotiles = unique([hitstruct.tile_id]');

for i=1:length(dotiles)
	
% Construct the path to the tile entity data

% 	tileid = find(dotiles(i)==[tFT.fac_id]);
	tileid = dotiles(i); % similar change to others?	

	hits = find([hitstruct.tile_id] == dotiles(i));
	
	tilepath = tFT(tileid).tile_name;
	tilepath = [ strrep(tilepath,'\',filesep) filesep];
	if isunix; tilepath = lower(tilepath); end
	tilepath = [themepath tilepath];
	tilepath = strrep(tilepath,[filesep filesep],filesep); % for the Browse library

% Read the face table for that tile
				
	EDG = vmap0read(tilepath,'EDG');
	FAC = vmap0read(tilepath,'FAC');
	RNG =  vmap0read(tilepath,'RNG');

% Construct the patch for the hits on that tile
	
	for n=1:length(hits)
		[lat,lon] = vmap0gbfac(hitstruct(hits(n)).fac_id,EDG,FAC,RNG);
		[latm,lonm] = meanm(lat(~isnan(lat)),lon(~isnan(lon)));
		
		[description,VDT,cVDT] = descript(hits(n),hitstruct,FTfield,FTfilename,themepath,VDT,cVDT);
		
		m=m+1;
		savestruct(m).string = hitstruct(hits(n)).nam;
		savestruct(m).lat = lat;
		savestruct(m).long = lon;
		savestruct(m).meanlat = latm;
		savestruct(m).meanlong = lonm;
		savestruct(m).type = 'patch';
		savestruct(m).tag = description;
		savestruct(m).altitude = [];
		savestruct(m).otherproperty = {};
	end

end % for i
