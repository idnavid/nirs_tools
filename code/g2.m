function [v a] = g2(t)
%UNTITLED2 Summary of this function goes here
%   t: Argument for g1
%   v: Value of g1 at t

v = zeros(size(t));
v = t.^12.*exp(-t./0.7).*(t>=0);
a = max(v);
%--------------------
%b1 = 0;
%b2 = 0;
%c = 0.4;

%v = c.*((t-b1).^12).*exp(-(t-b2)./0.7);
%a = max(v);

end

