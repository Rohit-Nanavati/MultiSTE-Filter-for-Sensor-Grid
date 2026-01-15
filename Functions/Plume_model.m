function [ C ] = Plume_model(psource,p,truesource)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    %     global s

    phi = truesource.phi;
    u = truesource.u;
    D = truesource.ci;
    t = truesource.cii;

    if isstruct(p) == 1
        px = p.x_matrix;
        py = p.y_matrix;
        C = zeros(size(px));
    else
        px = p(:,1);
        py = p(:,2);
        C = zeros(size(px));
    end
    % Loop over sources represented in the particle state, not truth.
    % This avoids dimension mismatch when the PF uses a different
    % source count than the ground truth (common in RJMCMC/unknown-K).
    N_sources = size(psource.x, 2);
    lamda = sqrt((D.*t)./(1+ (u.^2.*t)./(4*D)));

    for ns = 1:N_sources
        sx = psource.x(:,ns);
        sy = psource.y(:,ns);
        Q  = psource.Q(:,ns);

        module_dist = sqrt((sx-px).^2 + (sy-py).^2)+0.01;

        windcomponent = (px - sx).*cos(phi) + (py - sy).*sin(phi);

        C = C + Q./(2*pi.*D).*exp((windcomponent.*u)./(2.*D)).*besselk(0,module_dist./lamda);
    end
        C(besselk(0,module_dist./lamda)<=1e-30) = 0;
end

