%Calculates discriminability index " d' " and criterion "c" according to animal's Hit (H),
%Misses (M), False alarms (FA) and Correct rejects (CR) total numbers
%by S.M. 11/1/2019
%Updated by S.M. 12/03/2020
function [dp, c]=dprime_criterion_calc(a,b,c,d)

if nargin==4
    H=a;
    M=b;
    CR=c;
    FA=d;
    
    %calculate H and FA probabolity from values
    total_go=H+M;
    total_nogo=FA+CR;
    pH=H/total_go;
    pFA=FA/total_nogo;
    
elseif nargin==1
    H=a(1);
    M=a(2);
    CR=a(3);
    FA=a(4);
    
    %calculate H and FA probabolity from values
    total_go=H+M;
    total_nogo=FA+CR;
    pH=H/total_go;
    pFA=FA/total_nogo;
    
elseif nargin==2
    
    pH=a;
    pFA=b;
    total_go=25;
    total_nogo=25;
end



%if probability of H or FA is 0 or 1 aproximate according to Kuchibhotla et
%al 2019 for zero=1/2n and for one=1-1/2n

%pH
if pH==0
    pH=1/(2*total_go);
elseif pH==1
    pH=1-(1/(2*total_go));
else
    pH=pH;
end 

%pFA
if pFA==0
    pFA=1/(2*total_nogo);
elseif pFA==1
    pFA=1-(1/(2*total_nogo));
else
    pFA=pFA;
end 

%convert H and FA to z-score
% zH=-1*norminv(pH);
% zFA=-1*norminv(pFA); this was before, changed on 2020/12/03 to below
zH=1*norminv(pH);
zFA=1*norminv(pFA);


%obtain d' score
% dp=zFA-zH; this was before, changed on 2020/12/03 to below
dp=zH-zFA;

%obtain criterion score (z-score of FA)
c=-(zH+zFA)/2; %previously: c=zFA changed on 2020/12/03
end 
