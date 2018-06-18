function s = drift_mri(n)
% Return a quadratic and sinusoid drift function
%   n: length of the drift
%   p: polynomial drift
%   s: sinusoid drift

t = 1:n;
t = t/n;

%p = -7.8737+(47.5836.*t)-((32.6734.*t).*t);
a = 1;
s = a*sin(pi*(t-0.21));

end

