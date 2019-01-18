% The evaluation of associated Legendre transform
%
% The CURBF is credited to Eric Michielssen and Amir Boag, MULTILEVEL 
% EVALUATION OF ELECTROMAGNETIC FIELDS FOR THE RAPID SOLUTION OF SCAlTERlNG
% PROBLEMS, Microwave Opt. Technol. Lett., vol. 7, pp. 790-795, Dec. 1994.
%
% Copyright 2017 by Yingzhou Li and Haizhao Yang


function test_CUR_ALegendre ( )

tol = 1e-10;
NG = 10;  % number of Chebyshev pts
% set up problem size
N = 256;
kk = (0:(N-1))';
head = sprintf('root_%d.txt',N*2);
temp = load(head);
roots = temp(1:end/2).';
xx = roots.';
% set up order
order = 100;%120;
% ALegendrefun = @(x,k)fastALegendre(k,order,x);
    function mat = ALegendrefun(x,k)
        if isempty(x)
            mat = zeros(0,length(k));
            return;
        elseif isempty(k)
            mat = zeros(length(x),0);
            return;
        end
        mat = legendre(N,x,'norm');
        mat = mat(k+1,:).';
    end
ALegendremask = @(x,k)repmat(x,1,length(k)) ...
    > real(2*asin( sqrt(order^2-1/4)./(repmat(k.',length(x),1)+1/2) )/pi);

% Set up parameters
f = randn(N,1);
tic;
mFactor = CURBF_mask(@(x,k)ALegendrefun(x,k),xx,kk,ALegendremask,NG,tol);
mFactorT = toc;

tic;
yy = BF_mask_apply(mFactor,f);
ApplyT = toc;
RunT = mFactorT + ApplyT;

tic;
ytrue = ALegendrefun(xx,kk)*f;
relerr = norm(ytrue-yy)/norm(ytrue);
Td = toc;

% Plot patches

%subplot(1,2,1);
hold on;
for it = 1:length(mFactor)
    rectangle('Position', ...
        [min(kk(mFactor{it}.kidx)) min(xx(mFactor{it}.xidx))  ...
        max(kk(mFactor{it}.kidx))-min(kk(mFactor{it}.kidx)) ...
        max(xx(mFactor{it}.xidx))-min(xx(mFactor{it}.xidx))]);
end
axis tight;axis square;
xlabel('degree');ylabel('t');axis([0 N-1 0 1]);

curve = zeros(1,N);
for degree = 0:N-1
    curve(degree+1) = real(2*asin( sqrt(order^2-1/4)/(degree+1/2) )/pi);
end
%subplot(1,2,2);
plot(1:N,curve);axis tight;axis square;
xlabel('degree');ylabel('t');axis([0 N-1 0 1]);


disp(['------------------------------------------']);
disp(['N                 : ' num2str(N)]);
disp(['Chebyshev pts     : ' num2str(NG)]);
disp(['Tolerance         : ' num2str(tol)]);
disp(['Relative Error_2  : ' num2str(relerr)]);
disp(['Direct Time       : ' num2str(Td) ' s']);
disp(['Running Time      : ' num2str(RunT/60) ' mins']);
disp(['Factorization Time: ' num2str(mFactorT/60) ' mins']);
disp(['Applying Time     : ' num2str(ApplyT) ' s']);
disp(['------------------------------------------']);

end