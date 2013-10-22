function strout1=bldstr1(ordbig,Jbig,Dbig,i);
% strout1: subfunction of surrgt2.m; build ascii string summarizing surrogate tree info
% CALL: strout1=bldstr1(ordbig,Jbig,Dbig,i);
%
% Meko 4-21-97
%
%******************* IN ***** (see calling function for details on input) ***
%
% ordbig (1 x ?)i order of each time series
% Jthis (? x 2)i  like jthis, but for all series
% Dthis (? x 3)i  like dthis, but for all series
% i (1 x 1)i row index of Jthis and Dthis for current key series
%
%***************** OUT **********************
%
% strout1 (1 x ?)c   summary of daisy chain information; one line; several possible forms:
%
% 6(0.67/235)-12    Series 6 daisy chained to zero-order series 12. r=0.67 based on 235 yr.
% 5(O.43/54)-23-11  series 5(r=0.43,n=54yr) sieth series 23. Series 23 in turn daisy
%   chained to series 11
% 8*(0.12/32)-56  series 8 assigned the surrogate series 56. Correlation for the 32 yr
%   of overlap 0.12
% 34*(NaN/NaN)-2  series 34 assigned the surrogate series 2. No overlap of minimum
%   number of years (mlap in calling function) for correlation computation

dthis=Dbig(i,:);
jthis=Jbig(i,:);
order=ordbig(jthis(1));

% Set up a few formats
fmt1='%3.0f';
fmt2='%3.0f*';
fmt3='(%5.2f/%4.0f)-';
fmt4='%3.0f-';
fmt5='%3.0f\n';

strord=int2str(order);

% Take action depending on order of key series
switch strord
case '1'
   str1=sprintf(fmt1,jthis(1));
   str2=sprintf(fmt3,dthis(2),dthis(3));
   str3=sprintf(fmt5,jthis(2));
   strout1=[str1 str2 str3];
   return
case '99'
   str1=sprintf(fmt2,jthis(1));
   str2=sprintf(fmt3,dthis(2),dthis(3));
   str3=sprintf(fmt5,jthis(2));
   strout1=[str1 str2 str3];
   return
otherwise
   k=order;
   while k>0;
      if k==order;
         str1=sprintf(fmt1,jthis(1));
         str2=sprintf(fmt3,dthis(2),dthis(3));
         str3=sprintf(fmt4,jthis(2));
         strout1=[str1 str2 str3];
         
         % Find next link in daisy chain
         jrow=Jbig(:,1)==jthis(2);
         newj=Jbig(jrow,2);
         k=k-1;
      elseif k==1; % have reached last series in chain
         str1=sprintf(fmt5,newj);
         strout1=[strout1 str1];
         k=0;
         return;
      else
         str1=sprintf(fmt4,newj);
         strout1=[strout1 str1];
         jrow=Jbig(:,1)==newj;
         newj=Jbig(jrow,2);
         k=k-1;
      end; % if k==order
   end; % while k>0
end; % switch strord

      
  
