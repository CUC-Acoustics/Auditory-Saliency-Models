function r = LIN(in,L)
% This is the function corresponding to the lateral inhibition network, 
%
% INPUT
% -- in: input feature map
% -- L: size of the inhibitory field
%
% The coefficient distribution for lateral inhibition uses a hyperbolic distribution.

[xw,yw]=size(in);
smn=0;
for i=2:xw-1
    for j=2:yw-1
        for m = -L:L
            for n = -L:L
                d= (i-m)^2 + (j-n)^2;
                if d>=1
                    k= 0.7/d;
                else
                    k=0;
                end
                smn= smn + k*in(i+m,j+n);
            end
        end
        r(i,j)= in(i,j) - smn;
        smn=0;
    end
end
return
                
