clear
clc


for ii=1:3,
    fname = ['../test' num2str(ii) '.dat'];
    
    fid = fopen(fname, 'r');
    raw = fread(fid, inf, 'float32');
    fclose(fid);

    nIfg = (length(raw)-4)/2;
    wvl = raw(1);
    rng = raw(2);
    inc = raw(3) * pi / 180.0;
    delz = raw(4);

    bperp = raw(5:nIfg+4);
    ph =  raw(nIfg+5:2*nIfg+4);
    cph = exp(j * ph');

    maxK = 20.0/(wvl * rng * sin(inc)/4/pi);
    n_trial = (max(bperp) - min(bperp))*maxK/(2*pi);

    [K_r, C_r, coh_r] = ps_topofit(cph, bperp, n_trial, 'n');

    out = sprintf('file=%s K=%f C=%f coh=%f K_true=%f', fname, K_r,C_r, coh_r,delz/(wvl*rng*sin(inc)/4/pi));
    disp(out)

    figure('Name', num2str(ii));
    rph = C_r + K_r * bperp;

    rph = rph - round( rph/ (2*pi)) * 2*pi;

    scatter(bperp, ph, 'r');
    hold on;
    scatter(bperp, rph,'k');
end
