function com2html
% com2html:  convert comments of functions into html code as required by my web page
% com2html;
% Last revised 5-5-99
% 
% com2html is used to maintain my "tree-ring toolbox" web page.  com2html.m can be used
% optionally to cull the comments for a single function, or for all tree-ring toolbox
% functions at once.
%
%******** IN
%
% No args
% User prompted whether to handle one or many functions
% If a single function, user is promted to click on that m-file.  If many files, user
%   should make an ascii file named toolfls.txt with names of the functions, one per
%   line, in the directory c:\mlb\manip
%
%******** OUT
%
% No args
% User prompted for name of destination file (which is then later inserted into
% the web page).
%
%****** NOTES
%
% Assumes the following form of the .m function:
%
% line1:  "function..."
% line2:  "help1", consisting of function name, colon, and brief word descrp
% line3:  calling syntax
% line4:  blank
% line5-  paragraph(s) of opening detailed description, each ending with <P>
% lines ?  %*** IN  ... input section
% lines?   %*** OUT ... output section
% lines?   %*** UW Files Called ... user-written files called
% lines?   %*** TOOLBOXES
% line?    %*** NOTES .. Only the comments above this line are extracted for the
%    html file
%
% 
% A <BR> is inserted at the end of each line, except for those in the 
% group of lines beginning with line 5 and ending with the %*** IN line.  Thes
% lies are simply copied as is, which means the .html file has no BR breaks
% in this section.
% 
% Leading "%" characters are removed from lines
%
% User should beware of unintended "<...>" occurrences in the comment section


%----- Prompt for file list or just one m-file
kmen1=menu('Choose one',...
   'Process multiple m-functions using file list',...
   'Process just one function, pointing to file');

%--- Optionally read the file-name file
if kmen1==1;
   [file3,path3]=uigetfile('c:\mlb\manip\toolfls.txt','File with m-function file names');
   pf3=[path3 file3];
   flist= textread(pf3,'%s','delimiter','\n','whitespace','');
   nfun = size(flist,1); % number of function files in list
else
   nfun=1;
end;

[file2,path2]=uiputfile('c:\web\mflhtml.txt','Outfile to write converted comments to');
pf2=[path2 file2];
fid2=fopen(pf2,'w');


%----- Loop over m-functions
for n1 = 1:nfun;
   if kmen1==2;
      % Prompt for input of name of .m function
      [file1,path1]=uigetfile('*.m','Input file whose comments to be converted');
      pf1=[path1 file1];
   else; %kmen1==1, meaning multiple functions from list
      pf1 = char(flist(n1));
      disp(['Working on ' pf1]);
      % Extract file name from path\file
      file1=pf1;
      slash = findstr(pf1,'\');
      if ~isempty(slash);
         path1 = file1(1:max(slash));
         file1(1:max(slash))=[];
      end
   end
   
            
   % Read the function in as cell array of strings
   scell = textread(pf1,'%s','delimiter','\n','whitespace','');
   
   % Convert to character
   s = char(scell);
   
   % Drop first line, which should be "function" line
   s(1,:)=[];
   [ms,ns]=size(s);
   
   % Find the row marking start of the "NOTES" 
   s1=s;
   s1(:,1:5)=[];  % lop off left cols so that NOTES is left aligned in s1
   rowlast = strmatch('NOTES',s1); % get the row beginning with this text
   if isempty(rowlast);
      error('No line with NOTES left aligned at column 6');
   end
   clear s1;
   
   % Lop off all rows above the NOTES
   s=s(1:(rowlast-1),:);
   
   % Check that all rows retained begin with %
   col1 = s(:,1);  % extract col 1
   L2 = strcmp(col1,{'%'});
   if ~all(L2);
      error('All lines above NOTES must start with % ');
   end
   
   % Drop the first column, so no longer starts with '%'
   s(:,1)=[];
   [ms,ns]=size(s);
   
   %---- Find the row number of the "IN" line.  All lines starting with this
   % will have <BR> appended
   s2=s(:,1:6);
   i3 = find(strcmp(s2,{'*** IN'}))-1;  
   if isempty(i3);
      error('No line with *** IN');
   end
   % lines 1-5 should have <BR> appended
   % lines 6-i3 are the paragraph section
   % lines (i3+1):ms should have <BR> appended
   
   nfirst=4;
   nsecond=i3-4;
   nthird=ms-i3;
   
   sfirst= [s(1:4,:)];
   ssecond = [s(5:i3,:)];
   sthird= [s((i3+1):ms,:)];
   
   fnpref=strtok(file1,'.');
   
   dq='"'; % double quote, needed in netscape html
   txt4=['<A name=' dq  strtok(file1,'''.''')     dq  '><H3>' fnpref '</H3> </A>'];
   txt5 = '<A href="toolbox.html#fcnlist"> Back to Function List       </A>';
   
   %------- Write to txt file
   fprintf(fid2,'%s\n','<HR>');
   fprintf(fid2,'%s\n',txt4);
   
   fprintf(fid2,'%s\n','<PRE>');
   for n = 1:nfirst;
      srow = deblank(sfirst(n,:));
      fprintf(fid2,'%s\n',srow);
   end
   fprintf(fid2,'%s\n','</PRE>');
   
   for n = 1:nsecond;
      srow = deblank(ssecond(n,:));
      fprintf(fid2,'%s\n',srow);
   end
   
   fprintf(fid2,'%s\n','<PRE>');
   for n = 1:nthird;
      srow = deblank(sthird(n,:));
      fprintf(fid2,'%s\n',srow);
   end
   fprintf(fid2,'%s\n','</PRE>');
   
   fprintf(fid2,'%s\n','<BR>');
   fprintf(fid2,'%s\',txt5);
   
end;
fclose (fid2);




