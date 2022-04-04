function [k_i] = computeK(Phase)
k_i=ones(length(Phase),1);
for i = 1:size(k_i)
    if Phase(i)==0;
    k_i(i)=.60;
elseif Phase(i)==1;
    k_i(i)= 2.4;
 
elseif Phase(i)==2;
    k_i(i)= 1.8;
 
elseif Phase(i)==3;
    k_i(i)= 1.1;
    
elseif Phase(i)==5;
    k_i(i)=1.5;
 
else %Phase==6;
    k_i(i)= 1.9;
 
end
end