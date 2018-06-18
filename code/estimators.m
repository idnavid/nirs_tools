function h = estimators()
h.lse = @ls;
h.dbe= @dbe;
h.rob = @rob;
h.clean_up = @clean_up;
h.reg_gaussprior = @lse_gaussianprior;
h.regdbe = @regularized_dbe;
h.bi_noise = @bi_noise;
h.huber = @huber_est;
end

%% 
function n = bi_noise(dim,L,pi1,mu1,mu2,sigma1,sigma2)
% generate contaminated noise
n1 = mu1+sqrt(sigma1)*randn(dim,floor(L*pi1)); 
pi2 = 1 - pi1;
n2 = mu2+sqrt(sigma2)*randn(dim,floor(L*pi2));
n = zeros(dim,L);
n(:,1:length(n1)+length(n2)) = [n1, n2]; 
n=n(:,randperm(length(n)));
n = n(:);
end


function b_hat = ls(X,y)
% least square estimator of X*b = y.
% uses pseudo-inverse of X.
b_hat=X\y;
end

function b_hat = huber_est(X,y)
n = length(y);
y1 = y(1:n-1,1);
y2 = y(2:n,1);
y_dbe = y1-y2;

p = [1 -1 zeros(1,n-2)];
r = [p(1) zeros(1,n-1)];
m_dbe = toeplitz(p',r);
m_dbe =m_dbe(:,1:n-1);
X_dbe = m_dbe'*X;
b_hat = robustfit(X_dbe,y_dbe,'huber');
b_hat = b_hat(2:end);
end

function b_hat = regularized_dbe(X,y)
n = length(y);
y1 = y(1:n-1,1);
y2 = y(2:n,1);
y_dbe = y1-y2;

p = [1 -1 zeros(1,n-2)];
r = [p(1) zeros(1,n-1)];
m_dbe = toeplitz(p',r);
m_dbe =m_dbe(:,1:n-1);
X_dbe = m_dbe'*X;


%% Find optimal regularization parameter lambda
k = 0.001:0.1:10;
cv_MSE = zeros(size(k));
for i = 1:length(k)
    regf=@(X,y,Xtest)(Xtest*(inv(X'*X+k(i)*eye(size(X,2)))*X'*y));
    cv_MSE(i) = crossval('mse',X_dbe,y_dbe,'Predfun',regf,'kfold',10);
end
cv_MSE = medfilt1(cv_MSE,3);
desired_MSE = mean(cv_MSE); 
[~,idx] = min(cv_MSE-desired_MSE);
lambda = k(idx);
b_hat = inv(X_dbe'*X_dbe + lambda*eye(size(X_dbe,2)))*X_dbe'*y_dbe;
end

function [b_hat,cv_vec] = lse_gaussianprior(X,y,cv_vec)
% Estimate b using a Gaussian prior as proposed in 
% Ciuciu, Philippe, et al. "Unsupervised robust nonparametric estimation of
% the hemodynamic response function for any fMRI experiment." 
% IEEE TMI (2003): 1235-1251.
y = y(:);
[n,k] = size(X);
if n~=length(y)
    X = X';
    [n,k] = size(X);
end
n = length(y);
y1 = y(1:n-1,1);
y2 = y(2:n,1);
y_dbe = y1-y2;

p = [1 -1 zeros(1,n-2)];
r = [p(1) zeros(1,n-1)];
m_dbe = toeplitz(p',r);
m_dbe =m_dbe(:,1:n-1);
X_dbe = m_dbe'*X;
X = X_dbe; 
y = y_dbe;

left_lag = circshift(eye(k),1);
left_lag(1,k) = 0;
right_lag = left_lag';
D = -2*eye(k) + left_lag + right_lag;
R = D'*D;

%% Find optimal regularization parameter lambda
if false
    disp('\tau not specified. Performing cross validation ...')
    tau_vals = 0:20:2000;
    cv_MSE = zeros(size(tau_vals));
    for i = 1:length(tau_vals)
        regf=@(X,y,Xtest)(Xtest*(inv(X'*X+tau_vals(i)*R)*X'*y));
        cv_MSE(i) = crossval('mse',X_dbe,y_dbe,'Predfun',regf,'kfold',10);
    end
    cv_MSE = medfilt1(cv_MSE,3);
    if isempty(cv_vec)
        cv_vec = cv_MSE;
    else 
        cv_vec = cv_vec + cv_MSE;
    end
    desired_MSE = 0*mean(cv_MSE);
    [~,idx] = min(cv_MSE-desired_MSE);
    tau = tau_vals(idx);
    disp(tau)
    %figure(100)
    %plot(tau_vals, cv_MSE)
    %title(['\tau = ',num2str(tau)]),pause(2),hold off
    %close(100)
else
    tau = 200;
end
% Gaussian prior results in an invertable covariance on b, defined R = D'D.

b_hat=inv(X'*X + tau*R)*X'*y;
end


function b_hat = dbe(X,y)
n = length(y);
y1 = y(1:n-1,1);
y2 = y(2:n,1);
y_dbe = y1-y2;

p = [1 -1 zeros(1,n-2)];
r = [p(1) zeros(1,n-1)];
m_dbe = toeplitz(p',r);
m_dbe =m_dbe(:,1:n-1);

X_dbe = m_dbe'*X;

b_hat = pinv(X_dbe'*X_dbe)*X_dbe'*y_dbe;

end



function b_rob = rob(X,y,bo,alpha)
n = length(y);
y1 = y(1:n-1,1);
y2 = y(2:n,1);
y_dbe = y1-y2;

p = [1 -1 zeros(1,n-2)];
r = [p(1) zeros(1,n-1)];
m_dbe = toeplitz(p',r);
m_dbe =m_dbe(:,1:n-1);

X_dbe = m_dbe'*X;

for i = 1:10
    w = abs(y_dbe - X_dbe*bo);
    w = exp(-alpha*w.^2);
    w = w/sum(w);
    W = diag(w);
    bo = pinv(X_dbe'*W*X_dbe)*X_dbe'*W*y_dbe;
end
b_rob = bo;

end


function b = clean_up(b)
b = b - (b(1));
end