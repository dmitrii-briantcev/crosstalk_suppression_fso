%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Abderrahmen Trichilli
%
% Generation of the Laguerre-Gaussian beam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LG=lg(p,l,params)

%Function to calculate the LG field with the following initial conditions;
%x SLM screen size
%y SLM screen size
%p radial mode index
%l azimuthal mode index
%w beam radius
%n=2 parameter
%Initial position
%ph initial beam phase

w = params.w0;
n = 2;
z = 0;
Lambda = params.lambda;
ph = 0;

yy=params.y;
xx=params.x;

xx=xx;                                                                     % 8e-6 size of one pixel in m
yy=yy;

[XX, YY]=meshgrid(xx,yy);
r=sqrt(XX.^2+YY.^2);   
w=w;                                                                       % Beam radius normalization to mm

        k = (2.*pi)./Lambda;
        zr=(pi.*w.^2)./Lambda;
        R=z.*(1+(zr./z).^2);

        Io=sqrt((2.*factorial(p))./(pi.*factorial(p + abs(l))));
        v= 1./(w);
        v1=((sqrt(2).*r)./w).^(abs(l));
        v2=exp(-((r.^n)./(w.^n)));
        phi=atan2(YY,XX)+ pi/2;
        rho=2.*((r.^2)./(w.^2));                                           % Radius of the beam curvature
      
        L=laguerre(p,abs(l),rho);                                          % Generalized Laguerre polynomials
        phase=exp(1i.*l.*phi+ -1i.*ph);

if z==0
    im1=1;  
    im2=1;
else
    im1=exp(-1i.*k.*r.^2)./(2.*R);
    im2=exp(-1i.*(2.*p+abs(l)+1).*atan(z./zr));
end
                
        LG=Io.*v.*v1.*v2.*L.*phase.*im1.*im2;
        LG=LG/norm(LG, 'fro');
end
        