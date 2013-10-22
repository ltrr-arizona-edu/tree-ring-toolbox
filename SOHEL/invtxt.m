function h = invtxt(fh,x,y)
%
% USAGE : function h = invtxt(fh,x,y)
% Puts a text string at a given location (x,y) in a particular 
% figure window fh. 
% Returns the text object handle h.
%___________________________________________________________________

  figure(fh);
  v=axis;
  hold on;
  gtlbl=[num2str(x),' (',num2str(y),')'];
  vx=(x-v(1))/(v(2)-v(1))+0.02;
  vy=(y-v(3))/(v(4)-v(3));
  h=text('Position',[vx,vy],'String',gtlbl,'Fontsize',8,...
         'Units','Normalized');


% End of file