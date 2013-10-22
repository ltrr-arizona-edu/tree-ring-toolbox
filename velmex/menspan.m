function [C3,Lout] = menuspan(C1,C2,Lin,nallow,strmenu)
% menuspan:  menu for changing span associated with a series
% [C3,Lout] = menucell (C1,C2,Lin,,nallow,strmenu);
% Last revised 2-17-02
%
% Menu for changing span associated with a series. Called by rwlspecs.m to allow series-by-series
% tailoring of the output period for an rwl series
%
%*** INPUT
%
% C1{} cell array of names of series in the rwlset
% C2{} cell array of 1x2 integer vectors of specified span for output ([] if no restriction)
% Lin(? x 1)L  logical vector, same size as C1, with 1 for series currently restricted, 0 for not, on input
% Lout(?x1)L ditto, on output
% nallow (1 x 2)i  minimum and maximum allowable numbers of series on output to be restricted
% strmenu (1 x ?)s   Question starting menu (e.g., 'Toggle to Y to accept series')
%
%
%*** OUTPUT
%
% C3{} cell array of 1x2 integer vectors of specified span for output ([] if no restriction)
% Lout (? x1)L indicator whether on output series is restricted (1) or not (0) on output
%
%*** REFERENCED -- NONE
%
%*** UW FUNCTIONS CALLED --NONE
%*** TOOLBOXES NEEDED -- NONE
%*** NOTES


nser = length(C1);
D1=C1;
D1{nser+1}='Return to previous menu';
D2=C2;
D2{nser+1}='';
D0=cellstr([char(D1) char(D2)]);
D=D0;
Lcurr=Lin;

mess1=['Number of selected series must be in range ' int2str(nallow)];

kwh1 = 1;
while kwh1==1;
    kmen1=menu(strmenu,D);
    if kmen1==nser+1;
        npicked = sum(Lcurr);
        if npicked<nallow(1) | npicked>nallow(2);
            uiwait(msgbox(mess1,'Message','modal'));
        else;
            Lout=Lcurr;
            C3=C1(Lout);
            kwh1=0;
        end;
    else; % You have toggled a series
        if Lcurr(kmen1)==1; % was on, now will be off
            Lcurr(kmen1)=0;
            D2{kmen1}='-N';
        elseif Lcurr(kmen1)==0; 
            Lcurr(kmen1)=1;
            D2{kmen1}='-Y';
        end;
        % Revise menu
        D=cellstr([char(D1) char(D2)]); 
    end;
   
end; 
