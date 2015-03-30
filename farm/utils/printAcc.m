% print accuracy
function acc = printAcc(gr,thresh,name,obs)
% gr can be growth or growth ratio
if nargin<2
    thresh=0;
end
if nargin<3
    name='';
end
if nargin<4
    obs=0;
end

%treat NaN as no-growth prediction
nNan=sum(isnan(gr));

if obs==0
    nwrong=sum(gr>thresh);
    accName='specificity';
else
    nwrong=sum(gr<=thresh)+nNan;
    accName='sensitivity';
end
n=length(gr);

fprintf('name: %s\n', name);
fprintf('There are %d NaN values\n', nNan);
fprintf('%d wrong out of %d, %s %d\n', nwrong, n, accName, (n-nwrong)/n);
