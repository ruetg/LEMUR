function acc = getacc(ero,FD)
    acc = zeros(size(ero));
    acc=ero;
    FD.ix=(FD.I);
    FD.ixc = FD.R(FD.I);
    for i = length(FD.I):-1:1
        if FD.ixc(i)~=FD.ix(i)
            acc(FD.ixc(i)) = acc(FD.ixc(i)) + acc(FD.ix(i));
        end
    end
end

