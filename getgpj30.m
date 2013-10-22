% getgpj30.m  get pointer to gcn pcp stations (out of 300) with
% more than 30 good consecutive years of data
%
% Assume you have run seasgcn.m with
% J=(1:300)'
% kopt=[1 2]
% T=[1837 1995;1837 1995]
% X ... as usual , from gcnp.mat
%
% save result in seasaout.mat
% clear workspace
% load seasaout.mat


N=G(:,2)-G(:,1)+1;

J = find(N>30);

file1=uiputfile('*.dat','J_.dat the station-pointer vector');

eval(['save ',file1,' J -ascii'])

