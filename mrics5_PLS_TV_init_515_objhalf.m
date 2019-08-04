function [u, iter, snr] = mrics5_PLS_TV_init_515_objhalf(Pr0, Pr, R, f, alpha, beta, gamma, mu, alphal, betal, mul, MaxIter, tol, mask)
% mrics5_PLS_TV_init_515_objhalf reconstructs the EAPs from 5-dimensional undersampled 
% (k,q)-space diffusion MRI data (2D slice) using CS-based framework with parallel level 
% sets regularizations 

% This implementation takes data undersampled based on conventional DSI 
% sampling scheme (515 q locations), 6-neighborhood is used (2 in x, y and 
% z directions respectively).

% Input
%   Pr0: EAP initializations 
%   Pr: ground-truth EAPs
%   R: 6-dimensional sampling mask
%   f: data in the (k,q)-space
%   lambda, beta, mu: parameters for surfacelet transform, parallel level set and TV 
%                     regularization terms respectively 
%   lambdal, betal, mul: lagrangian multipliers corresponding to urfacelet transform, 
%                        parallel level set and TV regularization terms respectively 
%   MaxIter: maximum number of iterations
%   tol: stopping criteria (difference between consecutive terations)
%   mask: q-space mask 

% Output
%   u: reconstructed EAP (6-dimensional matrix)
%   iter: number of iterations performed
%   snr: history of snrs 

% Author: Jiaqi Sun
% August 2017, Last revision: 08-20-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    % parameters for numerical optimization tool
    addpath('./fminlbfgs_version2c');
    options = struct('Display', 'off','Algorithm','quasi-newton','GradObj','on');

    % load 6-neighborhood for conventional DSI sampling scheme 
    load('neighbors_515.mat');

    snr = zeros(MaxIter,1);
    [rows,cols,~,~,~] = size(f);
    dimr = 16;
    dimsurf = 32;
    scale = sqrt(rows*cols*dimr^3);
    numcoeff = 6.125*(dimsurf^3);
    
    volume = rand(dimsurf,dimsurf,dimsurf);
    D = surflet(volume);
    clear volume
        
    %% Reserve memory for the auxillary variables
    f0 = f;
    u = fftshift(Pr0);
    clear Pr0
    ushift = ifftshift(u);
    s = zeros(rows,cols,dimr,dimr,dimr);
    for i = 1:rows
        for j = 1:cols
            s(i,j,:,:,:) = fftshift(fftn(fftshift(squeeze(ushift(i,j,:,:,:)))));
        end
    end
    s = real(s);
    s = s.*mask;

    dw = Dw(u);
    dx = Dx(u);
    dy = Dy(u);
    
    reconstLattice = zeros(rows,cols,dimr,dimr,dimr);
    f3dtS = zeros(rows,cols,dimr,dimr,dimr);
    f3dU = zeros(rows,cols,dimr,dimr,dimr);
    bw = zeros(rows,cols,numcoeff);
    bsx = zeros(rows,cols,dimr,dimr,dimr);
    bsy = zeros(rows,cols,dimr,dimr,dimr);
    bx = zeros(rows,cols,dimr,dimr,dimr);
    by = zeros(rows,cols,dimr,dimr,dimr);
    sx = zeros(rows,cols,dimr,dimr,dimr);
    sy = zeros(rows,cols,dimr,dimr,dimr);

    w = shrink1(dw+bw,alpha/alphal);
    x = shrink1(dx+bx,mu/mul);
    y = shrink1(dy+by,mu/mul);

    %% Build Kernels
    ker = zeros(rows,cols,dimr,dimr,dimr);
    ker(1,1,1,1,1) = 4;
    ker(2,1,1,1,1) = -1;ker(rows,1,1,1,1) = -1;
    ker(1,2,1,1,1) = -1;ker(1,cols,1,1,1) = -1;
    
    uker = ((R).*R) + gamma + alphal + mul*fftn((ker));
    sker = gamma + betal*fftn((ker));
    
    du = 1;
    iter = 0;
    
    %%  Do the reconstruction
    while (du>tol || iter ==1) && iter < MaxIter 
%         tic;
        u0 = u;
        iter = iter+1;
        disp(['Iteration ', num2str(iter)]);
        
        %% update u
        for i = 1:rows
            for j = 1:cols
                    reconst = D'*(squeeze(w(i,j,:)-bw(i,j,:)));  
                    reconstLattice(i,j,:,:,:) = reconst(9:24,9:24,9:24);
                    f3dtS(i,j,:,:,:) = ifftn(ifftshift(squeeze(s(i,j,:,:,:))));
            end
        end
        urhs = ifftn((R).*f)*scale + gamma*f3dtS + alphal*reconstLattice + mul*(Dxt(x-bx)+Dyt(y-by));
        u = ifftn(fftn(urhs)./uker);
        u = real(u);
        ushift = ifftshift(u);
        %% update s
        for i = 1:rows
            for j = 1:cols
                f3dU(i,j,:,:,:) = fftshift(fftn(fftshift(squeeze(ushift(i,j,:,:,:)))));
            end
        end      
        srhs = gamma*f3dU + betal*(Dxt(sx-bsx)+Dyt(sy-bsy));
        s = ifftn(fftn(srhs)./sker);
        s = real(s);
        s = s.*mask;

        %% update sx, sy, sz and w
        dw = Dw(u);
        w = shrink1(dw+bw,alpha/alphal);
        
        dsx = Dx(s);
        dsy = Dy(s);
        
        dxyzuv = cat(6,dsx+bsx,dsy+bsy);
        dxyzuv = real(reshape(dxyzuv,[rows*cols dimr*dimr*dimr 2]));   
        dxyzuv_cen = dxyzuv(:,ind,:);  
        clear dxyzuv
        
        x0 = cat(6,sx,sy);
        x0 = real(reshape(x0,[rows*cols dimr*dimr*dimr 2]));  
        x0_cen = x0(:,ind,:);  
        x0_cen = reshape(x0_cen,[rows*cols 2*size(x0_cen,2)]);
        clear x0
        
%         tic;
%         [xyzuv_cen,~,~,~] = fminunc(@(x_cen)objfun_half(x_cen,dxyzuv_cen,nei,beta,betal),x0_cen,options);
%         [xyzuv_cen,~,~,~] = fmincon(@(x_cen)objfun_half(x_cen,dxyzuv_cen,nei,beta,betal),x0_cen,[],[],[],[],-inf(size(x0_cen)),[],[],options);
        [xyzuv_cen,~,~,~] = fminlbfgs(@(x_cen)objfun_half(x_cen,dxyzuv_cen,nei,beta,betal),x0_cen,options);
%         toc;
        
        
        xyzuv_cen = reshape(xyzuv_cen,[rows cols size(ind,1) 2]);
        sx = reshape(sx,[rows cols dimr*dimr*dimr]);
        sy = reshape(sy,[rows cols dimr*dimr*dimr]);
        sx(:,:,ind) = squeeze(xyzuv_cen(:,:,:,1));
        sy(:,:,ind) = squeeze(xyzuv_cen(:,:,:,2));
        sx = reshape(sx,[rows cols dimr dimr dimr]);
        sy = reshape(sy,[rows cols dimr dimr dimr]);
        clear x0_cen xyzuv_cen
        
        for i = 4:14
            for j = 4:14
                for k = 10:14
                    sx(:,:,i,j,k) = sx(:,:,18-i,18-j,18-k);
                    sy(:,:,i,j,k) = sy(:,:,18-i,18-j,18-k);
                end
            end
        end
        sx = sx.*mask;
        sy = sy.*mask;
        
        fprintf(1,'\n');
        
        dx = Dx(u);
        dy = Dy(u);
        x = shrink1(dx+bx,mu/mul);
        y = shrink1(dy+by,mu/mul);
        
        %% update bregman parameters bsx bsy bsz bw
        bsx = bsx+dsx-sx;
        bsy = bsy+dsy-sy;
        bw = bw+dw-w;

        bx = bx+dx-x;
        by = by+dy-y;
        %% check against stopping criteria 
        du2 = abs(u-u0);
        du2 = du2.*du2;
        du = sum(du2(:))/(scale^2);
        f = f+f0-R.*fftn(u)/scale;
        disp(['du = ', num2str(du)]);
        
        Pr_recon1 = ifftshift(u);
        Pr_recon1(Pr_recon1<0) = 0;
        d1 = (Pr-Pr_recon1).^2;
        ss = Pr.^2;
        NMSE1 = mean(d1(:))/mean(ss(:));
        snr1 = 10*log10(1/NMSE1);
        snr(iter) = snr1;
        fprintf(1,'SNR is %.4f\n',snr1);
        
%         % stop if the snr decreases in 2 consecutive iterations
%         if iter>2 && snr(iter)<snr(iter-1) && snr(iter-1)<snr(iter-2)
%             break;
%         end
        clear du2 Pr_recon1 ss d1 snr1
    end
      
return;
end

%% functions
function d = Dw(u)
dimsurf = 32;
[rows,cols,~,~,~] = size(u); 
d = zeros(rows,cols,6.125*dimsurf^3);
for i = 1:rows
    for j = 1:cols
        volume = zeros(dimsurf,dimsurf,dimsurf);
        volume(9:24,9:24,9:24) = squeeze(u(i,j,:,:,:));
        % original surfacelet transform on zero-padded P(R) volume
        D = surflet(volume);
%             d(i,j,:) = D*volume(:);
        d(i,j,:) = D*volume;

    end
end
return
end

function ws = shrink1(w, lambda)
    ws = abs(w)-lambda;
    ws = ws.*(ws>0);
    ws = ws.*sign(w);
return;
end

function d = Dx(u)
[rows,cols,dimrx,dimry,dimrz] = size(u); 
d = zeros(rows,cols,dimrx,dimry,dimrz);
d(:,2:cols,:,:,:) = u(:,2:cols,:,:,:)-u(:,1:cols-1,:,:,:);
d(:,1,:,:,:) = u(:,1,:,:,:)-u(:,cols,:,:,:);
return
end

function d = Dxt(u)
[rows,cols,dimrx,dimry,dimrz] = size(u); 
d = zeros(rows,cols,dimrx,dimry,dimrz);
d(:,1:cols-1,:,:,:) = u(:,1:cols-1,:,:,:)-u(:,2:cols,:,:,:);
d(:,cols,:,:,:) = u(:,cols,:,:,:)-u(:,1,:,:,:);
return
end

function d = Dy(u)
[rows,cols,dimrx,dimry,dimrz] = size(u); 
d = zeros(rows,cols,dimrx,dimry,dimrz);
d(2:rows,:,:,:,:) = u(2:rows,:,:,:,:)-u(1:rows-1,:,:,:,:);
d(1,:,:,:,:) = u(1,:,:,:,:)-u(rows,:,:,:,:);
return
end

function d = Dyt(u)
[rows,cols,dimrx,dimry,dimrz] = size(u); 
d = zeros(rows,cols,dimrx,dimry,dimrz);
d(1:rows-1,:,:,:,:) = u(1:rows-1,:,:,:,:)-u(2:rows,:,:,:,:);
d(rows,:,:,:,:) = u(rows,:,:,:,:)-u(1,:,:,:,:);
return
end

