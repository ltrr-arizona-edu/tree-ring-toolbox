function grplotd
% grplotd: plot groups of ring-width series together on figures
% CALL: gprlotd
% WHY: Usually run early in standardization of ring-width series.  Plots
%   help in decision of form of standarization curve to use. 
% WHO: Meko 4-14-99
%
%********************* IN ******************************************************
%
% No input arguments
%
% User prompted to point to files.  grplot1 can be run in two different modes.
%
% One mode is "single" for a single .mat ring-width file.  In this mode, user
% points to the desired .mat file, and optionally to a .mat file containing a vector
% of pointers to which series to plot and in what order.  Plots appear in figure
% widows with up to 8 pair of axes per window.  User can zoom on individual windows
% in this model.
%
% Another mode is "batch", in which the user points to a file that contains filenames
% of several .mat files of ring-width data (corresponding, say, to several sites).
%
% In either mode, series (cores) can be plotted in same order as they occur in the
% corresponding .mat file, or a pointer variable can be used to re-order plots or
% specify that some subset of all ring-files are to be plotted.
% The use of pointer variables is optional, in response to a prompt to the user.
% The pointer is a row-cell variable of matrices.  Each matrix goes with a 
% .mat file.  The matrices have as many rows as cores to be plotted from the .mat file,
% and 8 columns.  These columns are NaN-filled if needed.  
%
% Requirements for input files:
%  .mat ring-width file:  contains nms, yrs, X, as specified in rwlinp.m
%
%  .mat pointer file: contains a row-cell variable Vselect.  Each element of the 
%  cell is an n x 8 matric, where n is the number of figure windows to be produced.
%  The integer elements of the row are pointers to the series as ordered in the
%  string matrix nms in the .mat storage matrix.row of V tells which ring-width
%    to plot and in what order.  The col-size of 8 corresponds to the max number
%  of series that can be plotted in any figure window. For fewer than 8 plots,
%  NaN-fill the rows of the pointer matrix. Example for "single" mode:
%     V={[4 6 1 7 2 3 NaN NaN; 12 34 18 9 19 21 5 8]}  this specifies plotting
%     6 cores in figure window 1 and 8 in figure window 2.
%
% .txt file <batch mode only> with names of .mat files, one per line. Example:
%    shr.mat
%    dit.mat
%    ketnew.mat
%
%*************************** OUT *****************************************
%
% File output is optional, and consists of one or more postscript (.ps) files, one for
% each input .mat ring-width file.  The user can specify the filename in interactive 
% mode.  In batch mode, the .ps suffix is assigned to the same prefix as the input
% .mat file containing the ring-width file.  
%
% Interactive mode lets you view the figure windows as the program runs and afterwards.
% Batch mode does not save the figure windows, but overwrites them as it moves from
% ring-width file to ring-width file.
%
%****************************** UW FUNCTIONS CALLED ******************************
%
% eightplt.m - subfunction plots up to 8 sets of axes on page
%
%********************** NOTES *******************************************************8
%
%

% Close any open figure window
close all;

% Get current directory path
pathcurr=[eval('cd') '\']; 

%*************  What is mode: 'batch' or 'single'? *********

btchmode=questdlg('Batch mode (requires a .txt filelist)?');
switch btchmode;
case 'No';
   % no action needed
case 'Cancel';
   % no action needed;
case 'Yes';
   [file3,path3]=uigetfile('*.txt','Input .txt file with list of .mat filenames');
   pf3=[path3 file3];
   if ~(exist(pf3)==2);
      error(['Claimed file ' pf2 ' with list of .mat filenames does not exist']);
   end
   % Prompt for directory you want the .ps files to go into
   prompt={'Enter path:'};
   def={pathcurr};
   dlgTitle='Directory to store postcript (.ps) output files in';
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,def);
   pathps=answer{1};
   
   
end; % switch
% btchmode == 'Yes' means use batch mode
% pf3 would be the .txt filelist file if in batch mode


%*************** Compute number of ring-width files to treat in this run

if strcmp(btchmode,'No') | strcmp(btchmode,'Cancel'); % using 'single' mode
   nfiles=1; % number of .mat ring-width files
   filecell=[]; % row-cell holding names of the .mat files
else; % batch mode
   filecell = textread(pf3,'%s','delimiter','\n','whitespace','');
   nfiles=size(filecell,1); % number of .mat files to read
end; % strcmp(btchmode...)
% filecell is a row-cell with names of the .mat files
% nfiles is the number of .mat files


%****************** Plot all series in each file, or use a pointer to select files **

pointopt=questdlg('Use pointer to select ring-width series?');
switch pointopt;
case 'No'
   % no action needed;
case 'Cancel';
   % no action needed
case 'Yes';
   % load file with pointer variable;
   [file5,path5]=uigetfile('*.mat','Input file with selection pointer');
   pf5=[path5 file5];
   eval(['load ' pf5 ' Vselect;']); % pointer variable
   if exist('Vselect')~=1; error(['File ' pf5 ' does not contain Vselect']);end
   
   % Check size of Vselect
   [mcheck,ncheck]=size(Vselect);
   if (mcheck~=1 | ncheck~=nfiles);
      error('Vselect must be cell-row with row-size nfiles');
   end
   clear mcheck ncheck ;
   
end; % switch
% pointopt is 'Yes' if using a  pointer file
% If using a pointer file, Vselect is a rv or matrix of pointers to ring-width series

%------ Prompt for optional re-ordering of cores in .mat file from oldest to youngest for plotting 
oldfirst= questdlg('Arrange cores in order of age in plots?');
%  If pointopt=='yes', the ordering is within rows of Vmtx, which means
%  that retain the restriction that a specified group of cores is plotted
%  in each window
% If pointopt ~='yes', the ordering is over all cores in the .mat file
% Actual re-ordering is done within eightplt.m
 
%------ Loop over .mat files
for n = 1:nfiles;
   
   %-----  Get the .mat filename 
   if ~strcmp(btchmode,'Yes'); % single mode
      [flmat,path1]=uigetfile('*.mat','.MAT filename ?');
      pf1=[path1 flmat];
   else; % batch mode -- use a .txt list of .mat filenames
      pf1=filecell{n}; % file name, with or without drive\path
      len1=length(pf1);
      fslash=findstr(pf1,'\');
      if isempty(fslash); % no path prefix in file name; set path1 to current dir
         flmat=pf1; % filename (w/o path)
         path1=[eval('cd') '\']; % path
      else; % separate path and filename 
         flmat=pf1((max(fslash)+1):len1); % filename
         path1=pf1(1:max(fslash)); % path
      end
      % If needed, put on the .mat suffix
      if isempty(findstr(pf1,'.'));
         pf1=[pf1 '.mat'];
         flmat=[flmat '.mat'];
      end
      
   end
   % pf1 now is complete path\filename, including .mat
   % flmat  is the filename, with .mat
   % path1 is the path, with ending slash
   
      
   %----- Load the .mat file containing the ring-width series, names, years, etc;
   % Check that it holds the vital variables
   eval(['load ',pf1]);
   if ~all(exist('X')==1 & exist('nms')==1 & exist('yrs')==1);
      error('.mat input ring-width file does not have required variables');
   end
   ncore=size(nms,1);  %  number of ring-width series (cores) in the file,
   %  and maybe number to be plotted (depending on pointeropt)

   
   %--------  Make or pull matrix of selection pointers
   if strcmp(pointopt,'Yes'); % use read-in pointer variable
      Vmtx = Vselect{n}; % Matrix of pointers to core series
      Ltemp = ~isnan(Vmtx);
      ncore = sum(sum(Ltemp)); % total number of cores to be plotted
      numfig = size(Vmtx,1); % number of figure windows
      clear Ltemp;
   else; % not using pointer option;  build matrix of core pointers
      numfig = ceil(ncore/8); % number of figure windows
      i1=1:8;
      I1=repmat(i1,numfig,1);
      j1=(0:8:ncore-1)';
      J1=repmat(j1,1,8);
      Vmtx=I1+J1;
      ntemp = 8 * numfig;
      nquash = ntemp - ncore;
      if nquash>0;
         Vmtx(numfig,(8-nquash+1):8)=NaN;
      end
   end
   
      
   %------  Order cores from oldest to youngest before plotting.  
   if strcmp(oldfirst,'No'); % No re-ordering wanted
   else; % re-order
      if strcmp(pointopt,'Yes'); % using pointer option to specify cores in plot windows
         for ii=1:numfig; % loop over number of figure windows
            i1 = Vmtx(ii,:); % core sequence no's for this window
            L1 = isnan(i1); % any NaN's for this window
            if any(L1);
               nquash=sum(L1);
               i1(L1)=[];
               if length(i1)==1; % just one series to be plotted in window
                  Vmtx(ii,:)=[i1 repmat(NaN,1,7)];
               else; % more than one series; need to sort
                  yrborn = yrs(i1,1); % start year of rw for selected series
                  [yrtemp,isort]=sort(yrborn);
                  i2=i1(isort);
                  Vmtx(ii,:)=[i2,repmat(NaN,1,nquash)];
               end;
            else; % no NaN -- have 8 series to be plotted
               yrborn=yrs(i1,1); % start years
               [yrtemp,isort]=sort(yrborn);
               Vmtx(ii,:)=i1(isort);
            end; % if any(L1)
         end; % for ii=1:numfig
      else; % not using pointer option; re-order all cores in the .mat file
         i1=(1:ncore)'; % col vector of indices corresp to sequence of cores
         yrborn = yrs(:,1); % col vector of first years of ringwidths
         [yrtemp,isort]=sort(yrborn); 
         i2=i1(isort);  % core seq rearranged in order of start year of series
         if nquash>0; % need to make n of plot series multiple of 8
            i2=[i2; repmat(NaN,nquash,1)];
         end
         Vmtx=(reshape(i2,8,numfig))';
      end; % strcmp(pointopt,...);
   end; % strcmp(oldfirst,...) That's it for re-ordering
   
   
   
   % ---- Set datin{} contents that do not vary over figure windows
   datin{2}=X; % sov of ring widths
   datin{3}=nms; % core Ids
   datin{4}=yrs; % start year, end year, start row index of ring widths
   datin{7}=flmat; % (1 x ?)ch  name of .mat file from which ring widths come 
   %   (e.g., 'shr.mat') 
 
   %-------- Loop over figure windows
   for nn = 1:numfig;
      %-------  Load up input args for call to eightplt.m
      datin{1}=[nn numfig]; % window number & number of windows (not counting zoom)
      datin{5}=Vmtx(nn,:); % (? x 1)i row-index to nms of series to be plotted on this page
      datin{6}=[]; %  (1 x 2)d  specified first and last year for plot (as in zooming),
      %   or [], which means figure must cover 
      
      %---  Make the plots
      eightplt(datin);
   end; % for nn=1:numfig
   
   %------- ZOOM CAPABILITY
   
   if strcmp(btchmode,'Yes'); % zooming not implemented in batch mode
      % no zooming
   else; % 'single' mode
      kzoom=1;
      while kzoom==1;
         kmen1=menu('Choose One',...
            'Zoom a figure window',...
            'Quit zooming');
         if kmen1==1; % zoom a figure window
            
            if exist('nzw')==1;
               close(nzw);
            end
                          
            nzw = numfig+1; % figure window in which to put zoomed plot
            % close any existing zoom window
            
            % Build menu of possible windows to zoom
            figcell= cellstr(num2str((1:numfig)'));
            
            % Prompt for choosing a figure window , then switch to that window
            zwind=menu('Zoom on which window?',figcell);
            figure(zwind);
                       
            % Prompt user to point to zoom region
            hhlp1=helpdlg('Click L and R edge of desired zoom region');
            pause(1.5);
            xypnts=ginput(2);   % Get the corner points of the region to be zoomed in
            xmx=max(xypnts(:,1));
            xmn=min(xypnts(:,1));
                       
            % Make zoomed plot
            %-------  Load up input args for call to eightplt.m
            datin{1}=[nzw zwind]; % zoom window & source window
            datin{5}=Vmtx(zwind,:); % (? x 1)i row-index to nms of series to be plotted on this page
            datin{6}=[xmn xmx]; %  (1 x 2)d  specified first and last year for plot (as in zooming),
            %   or [], which means figure must cover 
            
            %---  Make the plots
            eightplt(datin);
            
         elseif kmen1==2; % quit zoom, but keep figures available
            kzoom=0;
         end; % if kmen1=1;
      end;  % while kzoom==1
   end; % strcmp(btchmode,...)
   
   %****************** OPTIONALLY MAKE POSTSCRIPT FILE OF THE FIGURE WINDOWS
   
   if ~strcmp(btchmode,'Yes'); % not in batch mode; optional making of .ps file
      psfile=questdlg('Make a postscript plot file?');
      switch psfile;
      case 'Yes'; %make a ps output file
         txt3=[' Output ps file for ' pf1];
         [file2,path2]=uiputfile('*.ps',txt3);
         pf2=[path2 file2];
         for npost = 1:numfig;
            figure(npost);
            eval(['print -dps -append pf2;']);
         end;
      case 'No';
      case 'Cancel';
      end; % switch psfile
      
   else; % batch mode;  assumed you want .ps files; these will go in current directory
      pf4=[pathps strtok(flmat,'.') '.ps'];
      for npost=1:numfig;
         figure(npost);
         eval(['print -dps -append ' pf4]);
      end
      close all 
   end; % strcmp(btchmode,...);
   
   % IF IN 'SINGLE' MODE, OPTIONALLY CLOSE ALL FIGURE WINDOWS
   if ~strcmp(btchmode,'Yes');
      kclose=questdlg('Close the figure windows?');
      switch kclose;
      case 'Yes';
         close all;
      case 'No';
      case 'Cancel';
      end
   end
         
end; % for n=1:nfiles
