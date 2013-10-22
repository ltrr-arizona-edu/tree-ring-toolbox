function u=splfres(p,fvec)
% u=splfres(p,fvec);
% Frequency response for spline of specified smoothing parameter
% Programmed by D. Meko   8-31-92
% Source:  Cook and Peters, 1981, Tree-ring Bulletin 41, 45-53
%
%*******   INPUT ARGUMENTS
%
% p  	smoothing parameter for spline
% fvec (3x1) starting frequency for freq response
%         ending   frequency 
%         number of frequency points resp desired at
%
%
%********    OUTPUT ARGUMENTS
%
% u (mux2)  freqencies (1/year) and corresp frequency response





%*********************  BEGIN CODE


f=linspace(fvec(1), fvec(2), fvec(3));   %  form vector of frequencies
w=2.0*pi*f;  % vector of angular frequencies

u=zeros(fvec(3),2); % allocate
u(:,1)=f';

costerm=cos(w');

term1=   p*(costerm+2);
term2=   6*(costerm-1) .^2;

% Remove any term2==0 to avoid division by 0
iout = find(term2==0);
if ~isempty(iout);
   nout=sum(iout);
   f(iout)=[];
   term2(iout)=[];
   term1(iout)=[];
   u(iout,:)=[];
end;


u(:,2) = 1 - (ones(length(f),1) ./   (1 + (term1 ./ term2)));
