
function temp = hamming_filter(temp)
temp=temp';
params.filterfunction=[.15 .7 .15];

MAXTP=size(temp,2);
runs_hammed=[];
a=1;
Run=[];
Run_hammed=[];
Run=temp(:,1+((a-1)*MAXTP):MAXTP*a);
for b = 1: size(temp,1)
    Run_hammed(b,:)=conv(Run(b,:),[params.filterfunction]/sum([params.filterfunction]));
end
if size([params.filterfunction],2) == 3
    Run_hammed=Run_hammed(:,2:(MAXTP+1));
    Run_hammed(:,1)=NaN; Run_hammed(:,MAXTP)=NaN;
elseif size([params.filterfunction],2 == 5)
    Run_hammed=Run_hammed(:,3:(MAXTP+2));
    Run_hammed(:,[1,2])=NaN; Run_hammed(:,[MAXTP-1, MAXTP])=NaN;
end
runs_hammed=[runs_hammed,Run_hammed];

temp=runs_hammed';
end