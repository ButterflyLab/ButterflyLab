clear; clc;
specFuncInit();

N = 2 .^ 13;%(9 : 15);
fun_m = @(n) round(1.0*n);

% Set up parameters
tol = 1e-10;
r = 128;
mR = 128;
nn = 256;

tstart = tic;
startT = char(datetime('now'));
cm = length(N);

for c = 1 : cm
    
    % set up problem size
    n = N(c);
    % set up order
    order = fun_m(n);
    isOdd = 1;
    if isOdd
        degreeVector = (1:2:(2*n-order))';
    else
        degreeVector = (2:2:(2*n-order))';
    end
    rootVector = ((n+1):(2*n))';
    
    %% visualize the transformation matrix, the discontinuity curve,
    %  and the matrix partition
    
    % compute quadrature points
    [rt,weight] = legenQuad(2*n);
    rt = acos(real(rt(:))); weight = real(weight(:));
    degree = (order:2*n-1)';
    % compute the transformation matrix using a fast algorithm
    % mat = legendre function(quadrature points, degrees) for
    % a fixed order when degree > order, mat = 0
    % in the nonzero region, odd columns are orthogonal to each
    % other even columns are orthogonal to each other
    ALegendrefun = @(x,k) ALegendrefun1(x,k,order,rt,degree,weight);
    
    mat = ALegendrefun(rootVector,degreeVector);
    
    % plot the complementary low-rank matrix
    pic = figure('visible', 'off');

    imagesc([order,2*n-1],[pi/2,0],mat);
    daspect([2 pi/2/n 1]);
    
    vec = degreeVector;
    % estimate and plot the discontiuity curve
    curve = zeros(1,2*n-order);
    for i = 1:2*n-order
        curve(i) = real(asin( sqrt(order^2-1/4)/(degree(i)+1/2) ));
    end
    hold on;
    plot(degree(vec),curve(vec),'r','LineWidth',2);
    axis([order 2*n-1 0 pi/2]); daspect([2 pi/2/n 1]);
    xlabel('degree'); ylabel('t');
    hold off;
    
    %% actual computation
    % apply the butterfly factorization to the complementary low-rank matrix
    
    % define function handle and avoid computing the whole matrix,
    % which is expensive
    fun = @(x,k) ALegendrefun(rootVector(x),degreeVector(k));
    
    % use the discontinuous curve to define the partition
    ALegendremask = @(x,k)repmat(rt(rootVector(x)),1,length(k)) ...
        > real(asin( sqrt(order^2-1/4)./(repmat(degree(degreeVector(k)).',...
        length(x),1)+1/2) ));
    
    proc = (c-1)/cm;
    fprintf(['%2.2f%% \t N = %i \t ', char(datetime('now','Format',...
        'MM-dd HH:mm:ss')),'\n'], proc*100, N(c));
        
    % apply hierarchical butterfly factorization
    mFactor = IDBF_mask(fun,rt(rootVector),degree(degreeVector),...
        ALegendremask,nn,r,mR,tol,'cheb',5,1,0);
    
    % visualize the matrix partitioning
    hold on;
    for it = 1:length(mFactor)
        rectangle('Position', ...
            [min(degree(degreeVector(mFactor{it}.kidx))) ...
            min(rt(rootVector(mFactor{it}.xidx))) ...
            max(degree(degreeVector(mFactor{it}.kidx))) - ...
            min(degree(degreeVector(mFactor{it}.kidx))) ...
            max(rt(rootVector(mFactor{it}.xidx))) - ...
            min(rt(rootVector(mFactor{it}.xidx)))]);
    end
    axis tight; axis square;
    axis([order 2*n-1 0 pi/2]); daspect([2 pi/2/n 1]);
    xlabel('degree'); ylabel('t');
    colorbar;
    hold off;
    saveas(gcf,['./results/',startT,' ',num2str(n),'.pdf']);
    close(pic);
    
end

T = toc(tstart) / 60;
fprintf('Finished. Running time = %.2f min\n', T);
