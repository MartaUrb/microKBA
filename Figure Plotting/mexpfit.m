function [mu,x0]=mexpfit(data,censor)
if isempty(censor)
    [cp,x] = ecdf(data);
else
    [cp,x] = ecdf(data,'censoring',censor);
end
x=(x(3:end)+x(2:end-1))/2;
cp=cp(2:end-1);
func = fittype('1-exp(-(x-x0)*k)');
f = fit(x,cp,func,'Startpoint',[1/nanmean(data) nanmin(data)]);
mu=1/f.k;
x0=f.x0;
