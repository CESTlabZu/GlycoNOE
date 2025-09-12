function  p = pulsesim1(dw, ksw1,ksw2,ksw3, ksw4,ksw5, ksw6, ksw7, kmw, mnots1,mnots2, mnots3,mnots4,mnots5, mnots6, mnots7, mnotw, mnotm,R1S, R2S1, R2S2, R2S3, R2S4, R2S5, R2S6, R2S7, R1W, R2W, R1M, R2M, sep1,sep2,sep3,sep4,sep5, sep6, sep7,duration, curve, angle, init)


    w1=getsatpulse(curve, angle, 1, duration);
    blah= pulsesolv3(w1, dw, ksw1, ksw2,ksw3,ksw4,ksw5, ksw6, ksw7, kmw, mnots1,mnots2,mnots3,mnots4,mnots5, mnots6, mnots7, mnotw, mnotm,R1S, R2S1, R2S2, R2S3, R2S4, R2S5, R2S6, R2S7,  R1W, R2W, R1M, R2M, sep1,sep2,sep3,sep4, sep5,sep6, sep7, init, duration);


   p = blah;%init(index(1), :);
end