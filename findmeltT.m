function [MeltT] = findmeltT(minP,maxP,minT,maxT,step,Pressure)
T1 = [minT:maxT]; %temperature range
P1 = [minP:maxP]; %pressure range
out = SF_WhichPhase({minP:step:maxP,minT:maxT});  
nT1=repmat(T1,(size(out,1)),1);
skip = 1;
for i = 1:size(out,1)
    idx(skip) = find((out(i,:)==0),1,'first');
    MeltT(skip) = nT1(skip,idx(skip));
    skip = skip+1;
end

end
% Melting_Temperature=MeltingTemps(1);
%  
% Melting_Temperature_list(i)= Melting_Temperature;

% 
% % [TT,PP] = meshgrid(T1,P1);
% nT1=repmat(T1,(size(out,1)),1);
% nP1=repmat(P1,(size(out,1)),1)';
% skip = 1;
% for i = 1:length(row)
%     meltT(skip) = nT1(row(i),col(i));
%     meltP(skip) = nP1(row(i),col(i));
%     skip = skip+1;
% end

% 
% [MeltT,MeltP] = (TT(x),PP(y));
% mm = nT1(y,x);

% mT=nT1(find(out==0));
%  
% P1_short=minP:step:maxP;
%  
% nP1=reshape(repmat(P1_short',(size(out,2)),1)', (size(out,1)),(size(out,2)));
%  
% mP=nP1(find(out==0));
% 
% tol = 5;
% skip = 1;
% for i = 1:length(Pressure)
% MeltingTemps=mT(abs(mP-Pressure(i)) < tol);
% Melting_Temperature(skip)=MeltingTemps;
% skip = skip+1;
% end
% MeltT= Melting_Temperature;