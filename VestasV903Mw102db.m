function [V,P,Prated] = VestasV903Mw102db(plotme);
P = [0
0
0
77
190
353
576
822
1056
1282
1500
1713
1920
2120
2291
2404
2459
2481
2490
2495
2500
2507
2516
2525
2532]/1000;

V = 1:length(P);
Prated = 2.5; %[Mw]
if nargin>0
    plot(V,P);
end
