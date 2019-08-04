function [f, g] = objfun_half(x,A,nei,beta,betal)
    sizex = [size(x,1) length(nei) 2];
    x = reshape(x,sizex);
    g = zeros(sizex);
    dxA = x-A;
    normX = sqrt(sum(x.*x,3));
    for i = 1:length(nei)
        nei1 = nei(i,1);
        nei2 = nei(i,2);
        nei3 = nei(i,3);
        x0 = squeeze(x(:,i,1));
        y0 = squeeze(x(:,i,2));
        a0 = squeeze(dxA(:,i,1));
        b0 = squeeze(dxA(:,i,2));
        norm0 = squeeze(normX(:,i))+1e-10;
        
        if nei1 > 0
            xx = squeeze(x(:,nei1,1));
            yx = squeeze(x(:,nei1,2));
            ax = squeeze(dxA(:,nei1,1));
            bx = squeeze(dxA(:,nei1,2));
            normx = squeeze(normX(:,nei1))+1e-10;
            f = sum(beta*(norm0.*normx-x0.*xx-y0.*yx)+betal/2*(a0.*a0+b0.*b0+ax.*ax+bx.*bx));
            g(:,i,1) = g(:,i,1) + beta*(normx./(norm0+1e-10).*x0-xx)+betal*a0;
            g(:,i,2) = g(:,i,2) + beta*(normx./(norm0+1e-10).*y0-yx)+betal*b0;
            g(:,nei1,1) = g(:,nei1,1) + beta*(norm0./(normx+1e-10).*xx-x0)+betal*ax;
            g(:,nei1,2) = g(:,nei1,2) + beta*(norm0./(normx+1e-10).*yx-y0)+betal*bx;
        end
       
        if nei2 > 0
            xy = squeeze(x(:,nei2,1));
            yy = squeeze(x(:,nei2,2));
            ay = squeeze(dxA(:,nei2,1));
            by = squeeze(dxA(:,nei2,2));
            normy = squeeze(normX(:,nei2))+1e-10;
            f = f + sum(beta*(norm0.*normy-x0.*xy-y0.*yy)+betal/2*(a0.*a0+b0.*b0+ay.*ay+by.*by));
            g(:,i,1) = g(:,i,1) + beta*(normy./(norm0+1e-10).*x0-xy)+betal*a0;
            g(:,i,2) = g(:,i,2) + beta*(normy./(norm0+1e-10).*y0-yy)+betal*b0;
            g(:,nei2,1) = g(:,nei2,1) + beta*(norm0./(normy+1e-10).*xy-x0)+betal*ay;
            g(:,nei2,2) = g(:,nei2,2) + beta*(norm0./(normy+1e-10).*yy-y0)+betal*by;
        end
        
        if nei3 > 0
            xz = squeeze(x(:,nei3,1));
            yz = squeeze(x(:,nei3,2));
            az = squeeze(dxA(:,nei3,1));
            bz = squeeze(dxA(:,nei3,2));
            normz = squeeze(normX(:,nei3))+1e-10;
            f = f + sum(beta*(norm0.*normz-x0.*xz-y0.*yz)+betal/2*(a0.*a0+b0.*b0+az.*az+bz.*bz));
            g(:,i,1) = g(:,i,1) + beta*(normz./(norm0+1e-10).*x0-xz)+betal*a0;
            g(:,i,2) = g(:,i,2) + beta*(normz./(norm0+1e-10).*y0-yz)+betal*b0;
            g(:,nei3,1) = g(:,nei3,1) + beta*(norm0./(normz+1e-10).*xz-x0)+betal*az;
            g(:,nei3,2) = g(:,nei3,2) + beta*(norm0./(normz+1e-10).*yz-y0)+betal*bz;
        end
    end
    g = reshape(g,[size(g,1) size(g,2)*2]);
end

