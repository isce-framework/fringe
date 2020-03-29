clear;

nwin = 24;
npad = 8;
alpha = 0.2;
beta = 0.5;
lp = 30.0;

cJ = 1i;

indata = zeros(nwin);
for ii=1:nwin,
    for jj=1:nwin,
        indata(ii,jj) = exp(cJ * (0.0314*(ii-1) + 0.02 * (jj-1)));
    end;
end;


ph_out = filtwin(indata, alpha, beta, nwin, npad, lp);
angs = angle(ph_out);


resid = angle(indata) - angs;
sum( exp(cJ*resid(:)))/ (24*24)
