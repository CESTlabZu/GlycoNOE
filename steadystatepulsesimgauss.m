function p = steadystatepulsesimgauss(dw, ksw1,ksw2, ksw3,ksw4, ksw5, ksw6, ksw7, kmw, mnots1,mnots2, mnots3,mnots4,mnots5, mnots6, mnots7, mnotw, mnotm, R1S, R2S1, R2S2, R2S3, R2S4, R2S5, R2S6, R2S7, R1W, R2W, R1M, R2M, sep1,sep2,sep3,sep4,sep5, sep6, sep7, duration, curve, angle,  time, spoil)
p = zeros(time,26);
init = [0;0;mnots1;0;0;mnotw;mnotm;0;0;mnots2; 0; 0; mnots3; 0; 0; mnots4; 0; 0; mnots5; 0; 0; mnots6; 0; 0; mnots7;1];
t = 1;

steadystate = 1;
while (steadystate == 1 && t < time)
   blah = pulsesim1(dw, ksw1,ksw2, ksw3,ksw4,ksw5, ksw6, ksw7, kmw, mnots1,mnots2,mnots3, mnots4, mnots5, mnots6, mnots7, mnotw, mnotm, R1S, R2S1, R2S2, R2S3, R2S4, R2S5, R2S6, R2S7,R1W, R2W, R1M, R2M, sep1,sep2,sep3,sep4, sep5,sep6, sep7, duration, curve, angle, init); 
   %blah;
   p(t,:) = blah;
       if(spoil == 1)
            blah(4) = 0;
            blah(5) = 0;
            blah(1) = 0;
            blah(2) = 0;
            blah(8) = 0;
            blah(9) = 0;
            blah(11) = 0;
            blah(12) = 0;
            blah(14) = 0;
            blah(15) = 0;
            blah(17) = 0;
            blah(18) = 0;
            blah(20) = 0;
            blah(21) = 0;
            blah(23) = 0;
            blah(24) = 0;
       end
   init = blah;
  %t
   t = t+1;

end
p(t:time, :) = [];

%p = init;
end