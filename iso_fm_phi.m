function [phi, J] = iso_fm_phi(lnkf, lnCf, lnaf, lnhf, lnks, lnCs, lnas, lnRth, lnsx, lnsy, f, Nx, dr)
% Nu Nv N n Nf
arguments
    lnkf  (1,1,:,1,1)
    lnCf  (1,1,:,1,1)
    lnaf  (1,1,:,1,1)
    lnhf  (1,1,:,1,1)
    lnks  (1,1,:,1,1)
    lnCs  (1,1,:,1,1)
    lnas  (1,1,:,1,1)
    lnRth (1,1,:,1,1)
    lnsx  (1,1)
    lnsy  (1,1)
    f     (1,1,1,1,:)
    Nx    (1,1)
    dr    (:,1)
end
Df = exp(lnkf-lnCf); % mm^2/s
Ds = exp(lnks-lnCs); % mm^2/s

Lthf = sqrt(Df ./ pi ./ f); % um
Lths = sqrt(Ds ./ pi ./ f); % um

x_max = max(max(dr)*10, max(Lthf, Lths));
% 
% x_max = [1.160157547808308e+03	8.699954270147923e+02	6.524045328641100e+02	4.892343813370420e+02	3.668740295709943e+02	2.751167103297527e+02	2.063084279668762e+02	200	200	200	200	200	200	200	200	200	200
%          1.143557271246048e+03	8.575469757474357e+02	6.430695113435620e+02	4.822341027547891e+02	3.616245611984488e+02	2.711801643951071e+02	2.033564349657142e+02	200	200	200	200	200	200	200	200	200	200
%          1.128927369323667e+03	8.465760970127407e+02	6.348425129091256e+02	4.760647242686181e+02	3.569981799964890e+02	2.677108679215929e+02	2.007548296297684e+02	200	200	200	200	200	200	200	200	200	200
%          1.103631981425568e+03	8.276072321051381e+02	6.206178709573048e+02	4.653977476391202e+02	3.489990759909439e+02	2.617123861479865e+02	200	200	200	200	200	200	200	200	200	200	200
%          1.084447343212171e+03	8.132207830007736e+02	6.098295560810858e+02	4.573076527849884e+02	3.429323607068801e+02	2.571629914867984e+02	200	200	200	200	200	200	200	200	200	200	200
%          1.064748943663876e+03	7.984490672463900e+02	5.987523319749685e+02	4.490009065723331e+02	3.367031798236156e+02	2.524917648135539e+02	200	200	200	200	200	200	200	200	200	200	200
%          1.040652183181413e+03	7.803790460969196e+02	5.852017277524658e+02	4.388393869329224e+02	3.290831150880035e+02	2.467775223935799e+02	200	200	200	200	200	200	200	200	200	200	200];
% x_max = reshape(x_max, [1,1,7,1,17]);

% x_max = [1.103731805694436e+03	8.276820897463128e+02	6.206740062689461e+02	4.654398431842591e+02	3.490306431964822e+02	2.617360582126245e+02	2.014885275089084e+02	200	200	200	200	200	200	200	200	200	200];
% x_max = reshape(x_max, [1,1,1,1,17]);

dx = x_max ./ floor(Nx/2);
du = 1 ./ (Nx * dx);
steps = -floor(Nx/2) : ceil(Nx/2) - 1;
x = steps(:) .* dx;
y = steps(:).' .* dx;
u = steps(:) .* du;
v = steps(:).' .* du;

if nargout > 1
    [T0hat{1:nargout(@iso_fm_g)}] = iso_fm_g(lnkf, lnCf, lnaf, lnhf, lnks, lnCs, lnas, lnRth, lnsx, lnsy, f, u, v);
    T0hat = cat(ndims(T0hat{1})+1,T0hat{:});
    T0tilde = fftshift(ifft2(ifftshift(T0hat)));
    
    [~,~,N,No,Nf,Ng] = size(T0tilde);

    x_interp = sqrt(dr(:).^2/2);
    
    X_interp = [x_interp,x_interp];
    
    N_interp = size(X_interp,1);
    
    T0tilde_interp = zeros(N_interp, N, No, Nf, Ng);
    
    xq = X_interp(:,1);
    yq = X_interp(:,2);
    
    for nf = 1:Nf
    for no = 1:No
    for n = 1:N
    
        xvec = x(:,1,n,no,nf);
        yvec = y(1,:,n,no,nf);
        V    = T0tilde(:,:,n,no,nf,:);
    
        F = griddedInterpolant({xvec, yvec}, zeros(length(xvec), length(yvec)), 'linear', 'none');
        F.Values = V;
        T0tilde_interp(:,n,no,nf,:) = F(xq, yq);
    
    end
    end
    end
    
    T0tilde = T0tilde_interp(:,:,:,:,1);
    
    dT0tilde_dm = T0tilde_interp(:,:,:,:,2:end);
    
    phi = angle(T0tilde);        % phase
    g = (real(T0tilde).*imag(dT0tilde_dm) ...
              - imag(T0tilde).*real(dT0tilde_dm)) ...
              ./ (abs(T0tilde).^2);
    lengths = [length(lnkf), length(lnCf), length(lnaf), length(lnhf), length(lnks), length(lnCs), length(lnas), length(lnRth)];
    J = zeros(numel(phi), sum(lengths));
    Jind = 1;
    for i = 1:length(lengths)
        J(:,Jind:Jind+lengths(i)-1) = make_jac_i(lengths(i), g(:,:,:,:,i));
        Jind = Jind + lengths(i);
    end
else
    T0hat = iso_fm(lnkf, lnCf, lnaf, lnhf, lnks, lnCs, lnas, lnRth, lnsx, lnsy, f, u, v);
    T0tilde = fftshift(ifft2(ifftshift(T0hat)));
    
    [~,~,N,No,Nf] = size(T0tilde);
    
    X_interp = [sqrt(dr(:).^2/2),sqrt(dr(:).^2/2)];
    
    N_interp = size(X_interp,1);
    
    T0tilde_interp = zeros(N_interp, N, No, Nf);
    
    xq = X_interp(:,1);
    yq = X_interp(:,2);
    
    for nf = 1:Nf
    for no = 1:No
    for n = 1:N
    
        xvec = x(:,1,n,no,nf);
        yvec = y(1,:,n,no,nf);
        V    = T0tilde(:,:,n,no,nf);
    
        F = griddedInterpolant({xvec, yvec}, zeros(length(xvec), length(yvec)), 'linear', 'none');
        F.Values = V;
        T0tilde_interp(:,n,no,nf) = F(xq, yq);
    
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
            gij(:,j,:,:) = gi(:,j,:,:);
            Ji(:,j) = gij(:);
        end
    end
end