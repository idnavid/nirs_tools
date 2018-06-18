function [v a] = g1(t)
%UNTITLED2 Summary of this function goes here
%   t: Argument for g1
%   v: Value of g1 at t

v = zeros(size(t));
v = t.^5.*exp(-t./0.9).*(t>=0);
a = max(v);
%------------------------
%b1 = 0;
%b2 = 0;

%v = ((t-b1).^5).*exp(-(t-b2)./0.9);

%a = max(v);
end

