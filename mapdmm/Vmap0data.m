function varargout = vmap0data(varargin)

%VMAP0DATA extracts selected data from the Vector Map Level 0 CD-ROMs
%
% struc = VMAP0DATA(library,latlim,lonlim,theme,topolevel) reads
% the data for the specified theme and topology level directly from the
% named VMAP0 CD-ROM. There are four CD's, one each for the library: 'NOAMER'
% (North America), 'SASAUS' (Southern Asia and Australia), 'EURNASIA'
% (Europe and Northern Asia) and 'SOAMAFR' (South America and Africa).
% The desired theme is specified by a two letter code string.  A list of
% valid codes will be displayed when an invalid code, such as '?', is
% entered.  The region of interest can be given as a point latitude and
% longitude or as a region with two-element vectors of latitude and
% longitude limits.  The units of latitude and longitude are degrees.
% The data covering the requested region is returned, but will include
% data extending to the edges of the 5 by 5 degree tiles.  The result is
% returned as a Mapping Toolbox Geographic Data Structure.
%
% [struc1, struc2, ...]  = ...
% VMAP0DATA(library,latlim,lonlim,theme,{topolevel1,topolevel2,...})
% reads several topology levels.  The levels must be specified as a cell
% array with the entries 'patch', 'line', 'point' or 'text'.  Entering
% {'all'} for the topology level argument is equivalent to {'patch',
% 'line', 'point','text'}. Upon output, the data are returned in the
% output arguments by topology level in the same order as they were
% requested.
%
% VMAP0DATA(devicename,library,latlim,...) specifies the logical device
% name of the CD-ROM for computers which do not automatically name the
% mounted disk.
%
% VMAP0 CD-ROMs are available from 
% 
% 	USGS Information Services (Map and Book Sales)
% 	Box 25286
% 	Denver Federal Center
% 	Denver, CO 80225
% 	Telephone: (303) 202-4700
% 	Fax: (303) 202-4693
%
% See also VMAP0GAZ, VMAP0READ, VMAP0RHEAD, MLAYERS, DISPLAYM, EXTRACTM

%  Copyright (c) 1996-97 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  W. Stumpf
%   $Revision: 1.6 $    $Date: 1997/03/28 23:11:49 $

% check for valid inputs

if nargin == 5
	library = varargin{1};
	latlim = varargin{2};
	lonlim = varargin{3};
	theme = varargin{4};
	topolevel = varargin{5};
	devicename = 'VMAP';
elseif nargin == 6
	devicename = varargin{1};
	library = varargin{2};
	latlim = varargin{3};
	lonlim = varargin{4};
	theme = varargin{5};
	topolevel = varargin{6};
else
	error('Incorrect number of input arguments')
end

%  Test the inputs

if  isequal(size(latlim),[1 1])
	latlim = latlim*[1 1];
elseif ~isequal(sort(size(latlim)),[1 2])
    error('Latitude limit input must be a scalar or 2 element vector')
end

if isequal(size(lonlim),[1 1])
	lonlim = lonlim*[1 1];
elseif ~isequal(sort(size(lonlim)),[1 2])
    error('Longitude limit input must be a scalar or 2 element vector')
end

%  Test for real inputs

if any([~isreal(latlim) ~isreal(lonlim)])
    warning('Imaginary parts of complex arguments ignored')
	latlim = real(latlim);   lonlim = real(lonlim);
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

for i=1:nargout
	varargout(i) = {[]};
end %for


% build the pathname so that [pathname filename] is the full filename

filepath = fullfile(devicename,'VMAPLV0',library,filesep);
if isunix; 
	filepath = fullfile(devicename,'vmaplv0',lower(library),filesep); 
end

% check that we can get into disk

dirstruc = dir(filepath);
if length(dirstruc) == 0; error(['Disk ' library ' not mounted']); end

% check for valid theme request

CAT = vmap0read(filepath,'CAT'); % was vmap0read

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

	tFT = vmap0read(filepath,'TILEREF.AFT'); % was vmap0read
	FBR = vmap0read(filepath,'FBR'); % was vmap0read

% find which tiles are fully or partially covered by the desired region

	dotiles = vmap0do(FBR,latlim,lonlim)	;

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

% read the Integer Value Description Table. This contains the integer
% values used to distinguish between types of features using integer keys

VDT = [];
if any(strcmp(uppercell(names),'INT.VDT'))==1
	VDT = vmap0read(themepath,'INT.VDT'); % was vmap0read
end


% read the Character Value Description Table if present.
% This contains the character values used to distinguish
% between types of features using character keys

cVDT = [];
if any(strcmp(uppercell(names),'CHAR.VDT'))==1
	cVDT = vmap0read(themepath,'CHAR.VDT'); % was vmap0read
end


% loop over points, lines, text
% handle faces separately

FACstruct = [];
EDGstruct = [];
ENDstruct = [];
TXTstruct = [];

for i=1:length(dotiles)

	FAC = [];
	EDG = [];
	END = [];
	TXT = [];
	SYM = [];


% extract pathname
% replace directory separator with the one for current platform

	tileid = find(dotiles(i)==[tFT.fac_id]);

	tilepath = tFT(tileid).tile_name;
	tilepath = [ strrep(tilepath,'\',filesep) filesep];
	if isunix; tilepath = lower(tilepath); end
	tilepath = [themepath tilepath];
	tilepath = strrep(tilepath,[filesep filesep],filesep); % for the Browse library

	dirstruc = dir(tilepath);
	tilenames = {dirstruc.name};
	tilenames = uppercell(tilenames);

%
% loop over feature topology level
%

	for j=1:4

%
% construct the feature table name. We know it will be of form
% '**'POINT.PFT', etc.
%
		switch j
% Patch
		case 1

			if any(topoindx ==1)

				EntityName =  'FAC';
%                if isunix; EntityName = 'FAC.'; end

				wildcard = '*.AFT';
				if strcmp(computer,'PCWIN')| isunix; wildcard = lower(wildcard); end

				FTfilenames = dir([themepath wildcard]);
				
				
				if isempty(FTfilenames) |  ~(any(strcmp(tilenames,EntityName))==1)
% 					warning('No patch data for this theme')
				else
				
					for k=1:length(FTfilenames)
	
						FTfilename = FTfilenames(k).name;
		
						if isempty(EDG); 
							EDGname =  'EDG';
%			                if isunix; EDGname = 'EDG.'; end
							EDG = vmap0read(tilepath,EDGname) ; % was ...
						end
	
						FACstruct = vmap0factl(themepath,tilepath,FTfilename,FACstruct,VDT,cVDT,EDG);
						if ~isempty(FACstruct);
							varargout(find(topoindx ==1))  = {FACstruct};
						end
	
					end % do k				
				
				end % if
					
			end % if
% Line		
		case 2

			if any(topoindx ==2)

				EntityName =  'EDG';
%                if isunix; EntityName = 'edg'; end

				wildcard = '*.LFT';
				if strcmp(computer,'PCWIN') | isunix; wildcard = lower(wildcard); end
				
				FTfilenames = dir([themepath wildcard]);
				if isempty(FTfilenames) |  ~(any(strcmp(tilenames,EntityName))==1)
% 					warning('No line data for this theme')
				else
				
					for k=1:length(FTfilenames)
	
						FTfilename = FTfilenames(k).name;

						if isempty(EDG); EDG = vmap0read(tilepath,EntityName); end; % was
	
						EDGstruct = vmap0edgtl(themepath,tilepath,FTfilename,EDGstruct,VDT,cVDT,EntityName,EDG);
						if ~isempty(EDGstruct);
							varargout(find(topoindx ==2)) = {EDGstruct};
						end

					end % do k				
				
				end % if
					
			end % if
			
% Point

		case 3

% 			if any(topoindx ==3)
% 
% 				FTfilename = [theme  'POINT.PFT'];
% 				EntityName =  'END';
% 	                        if isunix; EntityName = 'END.'; end
% 
% 				if ~(any(strcmp(names,FTfilename))==1)
% 					warning('No point data for this theme')
% 				elseif ~(any(strcmp(tilenames,EntityName))==1)
% 					warning('No point data for this theme')
% 				else
% 					if ~isempty(cVDT)
% 						cVDTindx = strmatch(FTfilename,strvcat(cVDT.table));
% 					end
% 					ENDstruct = vmap0endtl(themepath,tilepath,ENDstruct,VDT,FTfilename,'END',cVDT);
% 					if ~isempty(ENDstruct);
% 						varargout(find(topoindx ==3)) = {ENDstruct};
% 					end
% 				end %if
% 
% 			end % if

			if any(topoindx ==3)

				EntityName =  'END'; %new CND contains unique information?
%                if isunix; EntityName = 'END.'; end

				wildcard = '*.PFT';
				if strcmp(computer,'PCWIN') | isunix; wildcard = lower(wildcard); end
				
				FTfilenames = dir([themepath wildcard]);
				if isempty(FTfilenames) |  ~(any(strcmp(tilenames,EntityName))==1)
% 					warning('No point data for this theme')
				else
				
					for k=1:length(FTfilenames)
	
						FTfilename = FTfilenames(k).name;

						if isempty(END); END = vmap0read(tilepath,EntityName); end; % was
	
						EDGstruct = vmap0endtl(themepath,tilepath,FTfilename,EDGstruct,VDT,cVDT,EntityName,END);
						if ~isempty(EDGstruct);
							varargout(find(topoindx ==3)) = {EDGstruct};
						end

					end % do k				
				
				end % if
					
			end % if
% Text			
		otherwise

% 			if any(topoindx ==4)
% 				FTfilename = [theme 'TEXT.TFT'];
% 				EntityName =  'TXT';
%                                 if isunix; EntityName = 'TXT.'; end
% 
% 				if ~(any(strcmp(names,FTfilename))==1)
% 					warning('No text data for this theme')
% 				elseif ~(any(strcmp(tilenames,EntityName))==1)
% 					warning('No text data for this theme')
% 				else
% 
% 					TXTstruct = vmap0txttl(themepath,tilepath,TXTstruct,VDT,FTfilename,'TXT');
% 					if ~isempty(TXTstruct);
% 						varargout(find(topoindx ==4)) = {TXTstruct};
% 					end
% 				end
% 
% 			end % if

			if any(topoindx ==4)
				EntityName =  'TXT';
%                if isunix; EntityName = 'TXT.'; end

% check for and read SYMBOL.RAT file, which contains font sizes and colors

				SymbolName =  'SYMBOL.RAT';
				if strcmp(computer,'PCWIN') | isunix; SymbolName = lower(SymbolName); end
				SymbolDirListing = dir([themepath SymbolName]);
				if length(SymbolDirListing) ~= 0
					SYM = vmap0read(themepath,SymbolName); % was...
				end
					
				wildcard = '*.TFT';
				if strcmp(computer,'PCWIN') | isunix; wildcard = lower(wildcard); end
				
				FTfilenames = dir([themepath wildcard]);
				if isempty(FTfilenames) |  ~(any(strcmp(tilenames,EntityName))==1)
% 					warning('No text data for this theme')
				else
				
					for k=1:length(FTfilenames)
	
						FTfilename = FTfilenames(k).name;

						if isempty(TXT); TXT = vmap0read(tilepath,EntityName); end; % was
							
						TXTstruct = vmap0txttl(themepath,tilepath,TXTstruct,cVDT,FTfilename,TXT,SYM);
						if ~isempty(TXTstruct);
							varargout(find(topoindx ==4)) = {TXTstruct};
						end

					end % do k				
				
				end % if textdata
					
			end % if topoindx

		end % switch

	end % for j

end %for i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function carray = uppercell(carray)
%UPPERCELL converts a cellarray of strings to uppercase

for i=1:length(carray);carray{i} = upper(carray{i});end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function struc = vmap0edgtl(themepath,tilepath,FTfilename,struc,VDT,cVDT,EntityName,ET)

% Read the feature table and the entity table. The feature table contains
% the integer feature type. The entity file contains the coordinates.

etfields = fieldnames(ET); % match fieldnames to current entity table

% FTfilename = lower('DEPTHL.LFT')
% FTfilename = lower('polbndl.lft')
% FTfilename = lower('coastl.lft')

endofword = find(~isletter(FTfilename))-1;
matchthis = lower(FTfilename(1:endofword));
ftfieldindx = strmatch(matchthis,etfields); % which lft_id field goes with the current feature table filename

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Still need to make work for untiled layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isempty(ftfieldindx) % field name doesn't match what was expected for key
	indx = 1:length(ET); % must be in browse layer, so there are no tiles, and we want all elements
else
	eval(['indx = [ET.' etfields{ftfieldindx} '];']) % retrieve only elements in this tile
	ET = ET( indx > 0 );
	indx = indx(indx > 0); % These are the keys into the various feature tables. Positive ones go with the current feature table
	if isempty(indx); return; end % all of the LFT indices are keys into other LFT tables, so nothing to extract
end

[FT, FTfield]  = vmap0read(themepath,FTfilename, indx  ); % was

%
% find which fields of the feature table contain integer keys, text or values, for
% extraction as a tag. Later try to expand the text fields using the character value 
% description table (char.vdt). 
%

ftfields = fieldnames(FT);

FTindx = strmatch('int.vdt',{FTfield.VDTname}); 	% fields that are indices into the integer value description table


% Assume that 'crv' field is only occurance of a value. All other descriptive fields are 
% keys into integer or character description tables.

valfieldIDs = strmatch('crv',ftfields);

textfieldIDs = strmatch('T',{FTfield.type},'exact'); % into fields of FT

%
% Find out how many feature types we have by getting the indices
% into the feature table relating to the current topological entity
% (edge). Look for integer description keys, values, and text

valmat = [];
vdtfieldIDs = [];


if ~isempty(FTindx) % have integer feature, so gather  a matrix of the value combos

	valmat = [];
	
	for i=1:length(FTindx)
		
		eval(['valvec = [ FT.' ftfields{FTindx(i)} ']'' ;' ])
		valmat = [valmat valvec];
		
	end % for


end % if


% don't include VAL fields if they are keys into VDT?

valfieldIDs = valfieldIDs(find(~ismember( valfieldIDs, vdtfieldIDs)));

if  ~isempty(valfieldIDs) % have value fields, so gather them

	cells = struct2cell(FT);

	featvalvec = [cells{valfieldIDs,:}] ;
	featvalmat = reshape(featvalvec, length(valfieldIDs), length(featvalvec)/length(valfieldIDs));
	VALfeatvalmat = featvalmat';

	valmat = [valmat VALfeatvalmat];
end

if  ~isempty(textfieldIDs) % no integer feature, so depend on text strings in feature table. Only expect one.

	for i=1:length(textfieldIDs);
	 	eval(['strs = strvcat(FT.' ftfields{textfieldIDs(i)}	');']);
		valmat = [valmat double(strs)];
	end % for

end

[uniquecombos,cindx1,cindx2] = unique(valmat,'rows');
ncombos = size(cindx1,1);

%
% Extract the data
%

warned = 0;

i=length(struc);
for j=1:ncombos

	indx = find(cindx2==j);

	description = descript(valfieldIDs,textfieldIDs,ftfields,indx,FT,FTfield,FTindx,FTfilename,cVDT,VDT);

	struc(i+1).type= 'line';
	struc(i+1).otherproperty= {};
	struc(i+1).altitude= [];

	ll = [];
	for k=1:length(indx)
		llchunk = ET(indx(k)).coordinates;
		llchunk = llchunk(:,1:2); % for some reason we now have a trailing column of NaNs
		ll = [ll; NaN NaN; llchunk];
	end % for k


	struc(i+1).lat=ll(:,2);
	struc(i+1).long=ll(:,1);
	struc(i+1).tag=deblank(leadblnk(description));

	i=i+1;


end; %for j


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [struc,EDG] = vmap0factl(themepath,tilepath,FTfilename,struc,VDT,cVDT,EDG)

% read the feature table and the entity table. The feature table has in
% it the integer feature type and (optionally?) the name of the feature.
% The entity file contains the coordinates.

FAC = vmap0read(tilepath,'FAC'); % was
RNG = vmap0read(tilepath,'RNG');


etfields = fieldnames(FAC); % assume name of key into feature table is second field

endofword = find(~isletter(FTfilename))-1;
matchthis = lower(FTfilename(1:endofword));
ftfieldindx = strmatch(matchthis,etfields); % which lft_id field goes with the current feature table filename

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Still need to check for untiled layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isempty(ftfieldindx) % field name doesn't match what was expected for key
	indx = 1:length(FAC); % must be in browse layer, so there are no tiles, and we want all elements
else
	eval(['indx = [FAC.' etfields{ftfieldindx} '];']) % retrieve only elements in this tile
	FAC = FAC( indx > 0 );
	indx = indx(indx > 0); % These are the keys into the various feature tables. Positive ones go with the current feature table
	if isempty(indx); return; end % all of the FT indices are keys into other FT tables, so nothing to extract
end

[FT, FTfield]  = vmap0read(themepath,FTfilename,indx );
%
% find which fields of the feature table contain integer keys, text or values, for
% extraction as a tag. Later try to expand the text fields using the character value 
% description table (char.vdt). 
%

ftfields = fieldnames(FT);
FTindx = strmatch('int.vdt',{FTfield.VDTname}); 	% fields that are indices into the integer value description table


% Assume that 'crv' field is only occurance of a value. All other descriptive fields are 
% keys into integer or character description tables.

valfieldIDs = strmatch('crv',ftfields);
%>>>>>>>>>>>no crv elevation data for patches?<<<<<<<<<<<<<<

textfieldIDs = strmatch('T',{FTfield.type},'exact'); % into fields of FT


%
% Find out how many feature types we have by getting the indices
% into the feature table relating to the current topological entity
% (edge). Look for integer description keys, values, and text

valmat = [];
vdtfieldIDs = [];


if ~isempty(FTindx) % have integer feature, so gather  a matrix of the value combos

	valmat = [];
	
	for i=1:length(FTindx)
		
		eval(['valvec = [ FT.' ftfields{FTindx(i)} ']'' ;' ])
		valmat = [valmat valvec];
		
	end % for


end % if


% don't include VAL fields if they are keys into VDT?

valfieldIDs = valfieldIDs(find(~ismember( valfieldIDs, vdtfieldIDs)));

if  ~isempty(valfieldIDs) % have value fields, so gather them

	cells = struct2cell(FT);

	featvalvec = [cells{valfieldIDs,:}] ;
	featvalmat = reshape(featvalvec, length(valfieldIDs), length(featvalvec)/length(valfieldIDs));
	VALfeatvalmat = featvalmat';

	valmat = [valmat VALfeatvalmat];
end

if  ~isempty(textfieldIDs) % no integer feature, so depend on text strings in feature table. Only expect one.

	for i=1:length(textfieldIDs);
	 	eval(['strs = strvcat(FT.' ftfields{textfieldIDs(i)}	');']);
		valmat = [valmat double(strs)];
	end % for

end

[uniquecombos,cindx1,cindx2] = unique(valmat,'rows');
ncombos = size(cindx1,1);

%
% Extract the data
%

warned = 0;

i=length(struc);
for j=1:ncombos

	indx = find(cindx2==j);
	description = descript(valfieldIDs,textfieldIDs,ftfields,indx,FT,FTfield,FTindx,FTfilename,cVDT,VDT);

	for k=1:length(indx)
		
		[lat,long] = vmap0gbfac(FT(indx(k)).fac_id,EDG,FAC,RNG)   ;        % different
		
		if ~isempty(lat)

			struc(i+1).type= 'patch';                 % different
			struc(i+1).otherproperty= {};             % different

			struc(i+1).altitude= [];

			struc(i+1).lat=lat;
			struc(i+1).long=long;

			struc(i+1).tag=description;

			i=i+1;

		end %if
	end % for k

end; %for j

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function struc = vmap0endtl(themepath,tilepath,FTfilename,struc,VDT,cVDT,EntityName,ET)

% Read the feature table and the entity table. The feature table contains
% the integer feature type. The entity file contains the coordinates.

etfields = fieldnames(ET); % match fieldnames to current entity table

% FTfilename = lower('DEPTHL.LFT')
% FTfilename = lower('polbndl.lft')
% FTfilename = lower('coastl.lft')

endofword = find(~isletter(FTfilename))-1;
matchthis = lower(FTfilename(1:endofword));
ftfieldindx = strmatch(matchthis,etfields); % which lft_id field goes with the current feature table filename

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Still need to make work for untiled layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isempty(ftfieldindx) % field name doesn't match what was expected for key
	indx = 1:length(ET); % must be in browse layer, so there are no tiles, and we want all elements
else
	eval(['indx = [ET.' etfields{ftfieldindx} '];']) % retrieve only elements in this tile
	ET = ET( indx > 0 );
	indx = indx(indx > 0); % These are the keys into the various feature tables. Positive ones go with the current feature table
	if isempty(indx); return; end % all of the LFT indices are keys into other LFT tables, so nothing to extract
end

[FT, FTfield]  = vmap0read(themepath,FTfilename, indx  );

%
% find which fields of the feature table contain integer keys, text or values, for
% extraction as a tag. Later try to expand the text fields using the character value 
% description table (char.vdt). 
%

ftfields = fieldnames(FT);

FTindx = strmatch('int.vdt',{FTfield.VDTname}); 	% fields that are indices into the integer value description table


% Assume that 'crv' field is only occurance of a value. All other descriptive fields are 
% keys into integer or character description tables.

valfieldIDs = strmatch('crv',ftfields);

textfieldIDs = strmatch('T',{FTfield.type},'exact'); % into fields of FT

%
% Find out how many feature types we have by getting the indices
% into the feature table relating to the current topological entity
% (edge). Look for integer description keys, values, and text

valmat = [];
vdtfieldIDs = [];


if ~isempty(FTindx) % have integer feature, so gather  a matrix of the value combos

	valmat = [];
	
	for i=1:length(FTindx)
		
		eval(['valvec = [ FT.' ftfields{FTindx(i)} ']'' ;' ])
		valmat = [valmat valvec];
		
	end % for


end % if


% don't include VAL fields if they are keys into VDT?

valfieldIDs = valfieldIDs(find(~ismember( valfieldIDs, vdtfieldIDs)));

if  ~isempty(valfieldIDs) % have value fields, so gather them

	cells = struct2cell(FT);

	featvalvec = [cells{valfieldIDs,:}] ;
	featvalmat = reshape(featvalvec, length(valfieldIDs), length(featvalvec)/length(valfieldIDs));
	VALfeatvalmat = featvalmat';

	valmat = [valmat VALfeatvalmat];
end

if  ~isempty(textfieldIDs) % no integer feature, so depend on text strings in feature table. Only expect one.
	
	for i=1:length(textfieldIDs);
		
		eval(['indx = find(strcmp('''',{FT.' ftfields{textfieldIDs(i)} '}));' ] ) % empty strings
		if ~isempty(indx) & length(indx) < length(FT)
			for j=indx
				eval(['FT(j).' ftfields{textfieldIDs(i)} '= ''no string provided'';' ] ) % empty strings
			end
		end
		
		eval(['strs = strvcat(FT.' ftfields{textfieldIDs(i)}	');']);
		valmat = [valmat double(strs)];
	end % for

end

[uniquecombos,cindx1,cindx2] = unique(valmat,'rows');
ncombos = size(cindx1,1);

%
% Extract the data
%

warned = 0;

i=length(struc);
for j=1:ncombos

	indx = find(cindx2==j);

	description = descript(valfieldIDs,textfieldIDs,ftfields,indx,FT,FTfield,FTindx,FTfilename,cVDT,VDT);

	struc(i+1).type= 'line';
	struc(i+1).otherproperty= {};
	struc(i+1).altitude= [];

	ll = [];
	for k=1:length(indx)
		llchunk = ET(indx(k)).coordinate;
		llchunk = llchunk(:,1:2); % for some reason we now have a trailing column of NaNs
		ll = [ll; NaN NaN; llchunk];
	end % for k


	struc(i+1).lat=ll(:,2);
	struc(i+1).long=ll(:,1);
	struc(i+1).tag=deblank(leadblnk(description));

	i=i+1;


end; %for j


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function struc = vmap0txttl(themepath,tilepath,struc,cVDT,FTfilename,ET,SYM)


% read the feature table and the entity table. The feature table has in
% it the integer feature type and (optionally?) the name of the feature.
% The entity file contains the coordinates.

etfields = fieldnames(ET); % assume name of key into feature table is second field

endofword = find(~isletter(FTfilename))-1;
matchthis = lower(FTfilename(1:endofword));
ftfieldindx = strmatch(matchthis,etfields); % which lft_id field goes with the current feature table filename


if ftfieldindx % field name doesn't match what was expected for key
	indx = 1:length(ET); % must be in browse layer, so there are no tiles, and we want all elements
else
	eval(['indx = [ET.' etfields{ftfieldindx} '];']) % retrieve only elements in this tile
end


[FT,FTfield]  = vmap0read(themepath,FTfilename,indx);

%
% find out how many feature types we have by getting the indices
% into the feature table relating to the current topological entity
% (point, edge, face)
%


FTindx = strmatch(lower(FTfilename),{cVDT.table}); 	% the indices of the entries in the VDT pertinent to the current topology level

if isempty(FTindx)
	return % No appropriate Feature Table entry in Value Description Table
end

cVDT = cVDT(FTindx);
values = strvcat(cVDT.value);

% 
% vdtvalues = [VDT.value]	;					% all the values found in the Feature Table (over all topologies)
% featurevalues = unique(vdtvalues(FTindx));	% the values found in the Feature Table (for this topology level)
% 
% %
% % get the name of the attribute indexed in the Feature Table
% % from the Value Description Table, since this changes from
% % layer to layer
% %
% 
% attributeName = deblank(VDT(FTindx(1)).attribute);  % assume always same (questionable?)
% 
% %
% % Get all the attribute codes out of the Feature Table. Have to do
% % this by converting the structure to a cell array, to work around some
% % difficulty in 'eval'ing a big structure
% %
% 
% attributeCol = strmatch(deblank(attributeName), fieldnames(FT));
% cells = struct2cell(FT);
% allFeatVals = [cells{attributeCol,:}];

symcodes = [SYM.symbol_id];


i=length(struc);
for j=1:length(FT)

	struc(i+1).type= 'text';

	indx = strmatch(FT(j).f_code,values,'exact');
	
	struc(i+1).tag=deblank(cVDT(indx).description);
	struc(i+1).string=deblank(ET(j).string);
	struc(i+1).altitude= [];



	ll = ET(j).shape_line;


	struc(i+1).lat=ll(1,2);
	struc(i+1).long=ll(1,1);

	symindx = find(symcodes == FT(j).symbol_id);
	properties = {'fontsize',SYM(symindx).size};
	
	switch SYM(symindx).col
	case 1
		properties = { properties{:},'color','k'};
	case 4
		properties = { properties{:},'color','b'};
	case 9	
		properties = { properties{:},'color','r'};
	case 12
		properties = { properties{:},'color','m'};
	end	

	if length(ll(:,1)) > 1

		dx = ll(2,1) - ll(1,1);
		dy = ll(2,2) - ll(1,2);
		ang = 180/pi*(atan2(dy,dx));
		properties = { properties{:},'rotation',ang};

	end %if

	struc(i+1).otherproperty= properties;

	i=i+1;

end; %for j

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function description = descript(valfieldIDs,textfieldIDs,ftfields,indx,FT,FTfield,FTindx,FTfilename,cVDT,VDT)

% build up description codes

	description = '';

% extract integer based descriptions	

	for k = 1:length(FTindx)


		fieldname = ftfields{FTindx(k)};

		eval(['val = FT(indx(1)).' fieldname	';']);
		strng = num2str(val);
		
		tableindx = strmatch(lower(FTfilename), {VDT.table});
		fieldindx = strmatch(fieldname, {VDT.attribute});
		valueindx = find(val == [VDT.value]);
		VDTindx = intersect(tableindx, intersect(fieldindx,valueindx));
		
		if ~isempty(VDTindx)
			strng = VDT(VDTindx(1)).description;
			if strmatch(strng,'Unknown')  % expand labels of things with value unknown, 
				strng = [ FTfield(FTindx(k)).description ': ' strng];
			end
		end
	
		description = [description '; ' deblank(strng)];

	end %k
	
%  Extract text-based descriptions. In the VMAP0, it seems that all text is used 
%  as keys into the character description table. Look the text strings up in that
%  table, being careful to extract the correct occurance of duplicate codes. 
	
	if ~isempty(textfieldIDs) 

		for k = length(textfieldIDs):-1:1

			fieldname = ftfields{textfieldIDs(k)};

			eval(['strng = FT(indx(1)).' fieldname	';']);
			
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

		end % for k

	end %  ~isempty(textfieldIDs)

% extract value fields

	if ~isempty(valfieldIDs) % have contour values 

		for k = 1:length(valfieldIDs)

			eval(['val = FT(indx(1)).' ftfields{valfieldIDs(k)}	';']);
			description = [description '; ' num2str(val)];

		end %k

	end % ~isempty(valfieldIDs)


	description = description(3:length(description));
