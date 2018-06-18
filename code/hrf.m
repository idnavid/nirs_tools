function h = hrf(n,m)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

h = zeros(1,10*m);
k = round((9*(n^(1/5.5)-1))+1);

j = 1:k;

[h1 a1] = g1(1.5*(j-1)-5.5);
[h2 a2] = g2(1.5*(j-1)-5.5);


%a1 = max(h1);
%a2 = max(h2);

h = (h1/a1) - 0.4*(h2/a2); 

h(1,k+1:m) = 0;

%--------------------------------------
%t = 1:0.1:m;
%t = t/n;

%[v1 a1] = g1(1.5*(n*t-1));
%[v2 a2] = g2(1.5*(n*t-1));

%h = (v1./a1)-(v2./a2);

end


