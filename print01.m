function print01(S,wind1,wind2,nranked,fid2,yearend,monend);
% Dedicated subfunction of moving2.m to make tables of ranked results

%************  keep this while building subfunction
[file2,path2]=uiputfile('*.txt','File for summary of moving window analysis');
pf2=[path2 file2];
fid2=fopen(pf2,'w');
%***********************


%--------------- Cull desired window sizes

% Expand rv wind2 to matrix, same row size as wind1
W2=repmat(wind2,size(wind1,2),1);

% Make comparison matrix of wind1
W1=wind1'; % wind1 to cv
W1=repmat(W1,1,length(wind2));


% Logical pointer to desired windows
Ltemp=any((W1==W2)');  % a rv
wsizes=wind1(Ltemp); % desired window sizes

% Cull desired info on moving whatever, end year and end season/month
S1=S(1:nranked,Ltemp);
yearend=yearend(1:nranked,Ltemp);
monend=monend(1:nranked,Ltemp);

% Separate from previous printout
fprintf(fid2,'%s\n\n\n\n',blanks(5));

fmt1='%4.0f  %4.0f/%2.0f\t\t %g\n';
head2='Rank  End(Yr/mo)          Value';

%-------------  LOOP OVER SELECTED WINDOW SIZES
for n1 = 1:length(wsizes);
   width1=wsizes(n1);  % window width
   head1=sprintf('\n\n\n WINDOW SIZE  = %4.0f\n\n',width1);
   fprintf(fid2,'%s',head1);
   fprintf(fid2,'%s\n\n',head2);
   
   % Loop over ranks
   for n2=1:nranked;
      fprintf(fid2,fmt1,n2,yearend(n2,n1),monend(n2,n1),S1(n2,n1));
   end
  
      
end



fclose(fid2);