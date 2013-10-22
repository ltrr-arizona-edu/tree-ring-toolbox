function [yrs,strout]=respfn0h(yrsval)
% respfn0h:  subfunction of respfun1.m .  years for analysis


yrs=repmat(NaN,1,2);

txt1='Period for analysis (must be within ';
txt2=[int2str(yrsval(1)) '-' int2str(yrsval(2)) ')'];
txt3=[txt1 txt2];
txt4=['Start year earlier than ' int2str(yrsval(1)) ' or end year later than ' int2str(yrsval(2))];


kwh1=1;

while kwh1==1;
   prompt={'Enter start year','Enter ending year'};
   def = {int2str(yrsval(1)),int2str(yrsval(2))};
   titdlg = txt3;
   lineNo=1;
   answer=inputdlg(prompt,titdlg,lineNo,def);
   yrs(1)=str2num(answer{1});
   yrs(2)=str2num(answer{2});
   
   kwh1=0;
   
   if yrs(2)<=yrs(1);
      h1=helpdlg('Start year must precede end year','Hit OK and try again!');
      pause(3);
      if ishandle(h1);
         close(h1);
      end
      kwh1=1;
   end
   
   if yrs(1)<yrsval(1) | yrs(2)>yrsval(2);
      h1=helpdlg(txt4,'Hit OK and try that again');
      pause(3);
      if ishandle(h1);
         close(h1);
      end
      kwh1=1;
   end
   
end

strout='SELECTED PERIOD FOR ANALYSIS OF CLIMATE-TREE RING RELATIONSHIP';
strout=char(strout,...
   ['  ' int2str(yrs(1)) '-' int2str(yrs(2))]);   

