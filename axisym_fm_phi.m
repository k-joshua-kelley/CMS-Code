function [phi, J] = axisym_fm_phi(lnkf,lnCf,lnaf,lnhf,lnks_perp,lnks_par,lnCs,lnas,lnRth,vs1,vs2,vs3,lnsx,lnsy,f,Nx,X_probe,x_max)
% Nu Nv N n Nf
arguments
    lnkf      (1,1,:,1,1)
    lnCf      (1,1,:,1,1)
    lnaf      (1,1)
    lnhf      (1,1)
    lnks_perp (1,1,:,1,1)
    lnks_par  (1,1,:,1,1)
    lnCs      (1,1,:,1,1)
    lnas      (1,1)
    lnRth     (1,1)
    vs1       (1,1,1,:,1)
    vs2       (1,1,1,:,1)
    vs3       (1,1,1,:,1)
    lnsx      (1,1)
    lnsy      (1,1)
    f         (1,1,1,1,:)
    Nx        (1,1)
    X_probe   (:,2)
    x_max     (1,1,:,1,:)
end
dx = x_max ./ floor(Nx/2);
du = 1 ./ (Nx * dx);
steps = -floor(Nx/2) : ceil(Nx/2) - 1;
x = steps(:) .* dx;
y = steps(:).' .* dx;
u = steps(:) .* du;
v = steps(:).' .* du;

if nargout > 1
    [T0hat{1:nargout(@axisym_fm_g)}] = axisym_fm_g(lnkf,lnCf,lnaf,lnhf,lnks_perp,lnks_par,lnCs,lnas,lnRth,lnsx,lnsy,f,u,v,vs1,vs2,vs3);
    T0hat = cat(ndims(T0hat{1})+1,T0hat{:});
    T0tilde = fftshift(ifft2(ifftshift(T0hat)));
    
    [~,~,N,No,Nf,Ng] = size(T0tilde);
        
    Nprobe = size(X_probe,1);
    
    T0tilde_interp = zeros(N, No, Nf, Nprobe, Ng);
    
    xq = X_probe(:,1);
    yq = X_probe(:,2);
    
    for nf = 1:Nf
    for no = 1:No
    for n = 1:N
    
        xvec = x(:,1,n,1,nf);
        yvec = y(1,:,n,1,nf);
        V    = T0tilde(:,:,n,no,nf,:);
    
        F = griddedInterpolant({xvec, yvec}, zeros(length(xvec), length(yvec)), 'linear', 'none');
        F.Values = V;
        T0tilde_interp(n,no,nf,:,:) = F(xq, yq);
    
    end
    end
    end
    
    T0tilde = T0tilde_interp(:,:,:,:,1);
    
    dT0tilde_dm = T0tilde_interp(:,:,:,:,2:end);
    
    phi = angle(T0tilde);        % phase
    g = (real(T0tilde).*imag(dT0tilde_dm) ...
              - imag(T0tilde).*real(dT0tilde_dm)) ...
              ./ (abs(T0tilde).^2);
    lengths = [length(lnkf), length(lnCf), length(lnhf), length(lnks_perp), length(lnks_par), length(lnCs), length(vs1), length(vs2), length(vs3)];
    J = zeros(numel(phi), sum(lengths));
    Jind = 1;
    for i = 1:length(lengths)
        J(:,Jind:Jind+lengths(i)-1) = make_jac_i(lengths(i), g(:,:,:,:,i));
        Jind = Jind + lengths(i);
    end
else
    T0hat = axisym_fm(lnkf,lnCf,lnaf,lnhf,lnks_perp,lnks_par,lnCs,lnas,lnRth,lnsx,lnsy,f,u,v,vs1,vs2,vs3);
    T0tilde = fftshift(ifft2(ifftshift(T0hat)));
    
    [~,~,N,No,Nf] = size(T0tilde);
        
    Nprobe = size(X_probe,1);
    
    T0tilde_interp = zeros(N, No, Nf, Nprobe);
    
    xq = X_probe(:,1);
    yq = X_probe(:,2);
    
    for nf = 1:Nf
    for no = 1:No
    for n = 1:N
    
        xvec = x(:,1,n,1,nf);
        yvec = y(1,:,n,1,nf);
        V    = T0tilde(:,:,n,no,nf);
    
        F = griddedInterpolant({xvec, yvec}, zeros(length(xvec), length(yvec)), 'linear', 'none');
        F.Values = V;
        T0tilde_interp(n,no,nf,:) = F(xq, yq);
    
    end
    end
    end
    
    phi = angle(T0tilde_interp);
end
end

function Ji = make_jac_i(len, gi)
    Ji = zeros(numel(gi), len);
    if len == 1
        Ji = gi(:);
    else
        for j = 1:len
            gij = zeros(size(gi));
            dim = find(size(gi) == len);
            
            if isempty(dim)
                error('No matching dimension.');
            elseif numel(dim) > 1
                error('Multiple matching dimensions.');
            end
            
            switch dim
                case 1
                    gij(j,:,:,:) = gi(j,:,:,:);
                case 2
                    gij(:,j,:,:) = gi(:,j,:,:);
                case 3
                    gij(:,:,j,:) = gi(:,:,j,:);
                case 4
                    gij(:,:,:,j) = gi(:,:,:,j);
            end
            Ji(:,j) = gij(:);
        end
    end
end