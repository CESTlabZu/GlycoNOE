function l = runsteadysimgauss(ksw1,ksw2, ksw3,ksw4, ksw5,ksw6, ksw7, kmw, mnots1,mnots2, mnots3, mnots4, mnots5, mnots6, mnots7, mnotw, mnotm, R1S, R2S1, R2S2, R2S3, R2S4, R2S5, R2S6, R2S7, R1W, R2W, R1M, R2M, sep1,sep2,sep3,sep4, sep5, sep6, sep7, duration, curve, angle,time,  k, spoil)
l = zeros(89, 26);

for t=1:size(k,1)
    p = steadystatepulsesimgauss(k(t), ksw1,ksw2, ksw3,ksw4, ksw5, ksw6, ksw7, kmw, mnots1,mnots2, mnots3,mnots4, mnots5, mnots6, mnots7, mnotw, mnotm, R1S, R2S1, R2S2, R2S3, R2S4, R2S5, R2S6, R2S7, R1W, R2W, R1M, R2M, sep1,sep2, sep3,sep4,sep5, sep6, sep7, duration, curve, angle,  time, spoil);
    ind2 = size(p);
    l(t, :) = p(ind2(1), :);
end