function modes = BF_ModeDecom(img,opt)

VMDonly = opt.VMDonly;
N = size(img,1);
% some sample parameters for VMD
alpha2 = 2000;        % moderate bandwidth constraint
tau = 0;            % noise-tolerance (no strict fidelity enforcement)
DC = 0;             % no DC part imposed
init = 1;           % initialize omegas uniformly
tol = 1e-7;

modes = zeros([size(img), opt.numCom]);
for cntc = 1:size(img,2)
    fff = real(img(:,cntc));
    fff = fff(:).';
    [uc, u_hat, omega] = VMD(fff, alpha2, tau, opt.numCom, DC, init, tol,5);
    
    if ~VMDonly
        iniIF = zeros(opt.numCom,N);
        for k = 1:opt.numCom
            %% STFT
            window = 256;
            Nfrebin = 1024;
            [Spec,f] = STFT(uc(k,:)',N,Nfrebin,window);
            %% ridge extraction and smoothing
            bw = N/80;% the bandwidth of the TF filter for ridge extraction
            beta1 = 1e-4; % beta1 should be larger than the following beta
            delta = 20;
            [fidexmult, tfdv] = extridge_mult(uc(k,:), N, 1, delta, beta1,bw,Nfrebin,window);
            
            %% parameter setting
            alpha = 1e-5;
            beta = 1e-5; % this parameter can be smaller which will be helpful for the convergence, but it may cannot properly track fast varying IFs
            iniIF(k,:) = curvesmooth(f(fidexmult),beta);
        end
        var = 0;% noise variance
        [IFmset,IA,smset] = VNCMD(fff,N,iniIF,alpha,beta,var,tol);
        uc = smset(:,:,end).';
    end
    
    if ~VMDonly
        modes(:,cntc,:) = modes(:,cntc,:) + reshape(uc,[N,1,opt.numCom]);
    else
        modes(:,cntc,:) = modes(:,cntc,:) + reshape(uc.',[N,1,opt.numCom]);
    end
    
    if opt.isCmpl
        fff = imag(img(:,cntc));
        fff = fff(:).';
        [uc, u_hat, omega] = VMD(fff, alpha2, tau, opt.numCom, DC, init, tol,5);
        
        if ~VMDonly
            iniIF = zeros(opt.numCom,N);
            for k = 1:opt.numCom
                [Spec,f] = STFT(uc(k,:)',N,Nfrebin,window);
                [fidexmult, tfdv] = extridge_mult(uc(k,:), N, 1, delta, beta1,bw,Nfrebin,window);
                iniIF(k,:) = curvesmooth(f(fidexmult),beta);
            end
            [IFmset,IA,smset] = VNCMD(fff,N,iniIF,alpha,beta,var,tol);
            uc = smset(:,:,end).';
        end
        
        if ~VMDonly
            modes(:,cntc,:) = modes(:,cntc,:) + 1i*reshape(uc,[N,1,opt.numCom]);
        else
            modes(:,cntc,:) = modes(:,cntc,:) + 1i*reshape(uc.',[N,1,opt.numCom]);
        end
    end
end