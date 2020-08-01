thres=50;
[mm,nn]=size(Z1);
acc = zeros(mm,nn)+1;
FD.ix=fliplr(FD.I);
FD.ixc=FD.R(FD.ix);
if 1
[LON,LAT]=meshgrid(1:nn,1:mm);
end
for i=1:length(FD.ix)
    acc(FD.ixc(i))=acc(FD.ix(i))+acc(FD.ixc(i));
end

nplotx =zeros(1,nn*mm+1000);
nploty=zeros(1,nn*mm+1000);
lix=FD.ixc(i);
n=1;
for i=1:length(FD.ix)
    [y1,x1]=ind2sub([mm,nn],FD.ix(i));

    if (FD.ix(i)==lix&&acc(FD.ix(i))>thres)
        
        nploty(n)=LAT(FD.ix(i));
        
        nplotx(n)=LON(FD.ix(i));
    else
        nploty(n)=NaN;
        nplotx(n)=NaN;
        n=n+1;
        nploty(n)=LAT(FD.ix(i));
        
        nplotx(n)=LON(FD.ix(i));
    end
    n=n+1;
    lix=FD.ixc(i);
end
plot(nplotx,nploty,'k','linewidth',2);

