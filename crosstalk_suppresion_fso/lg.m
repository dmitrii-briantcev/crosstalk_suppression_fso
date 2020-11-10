%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Abderrahmen Trichilli
%
% Generation of the Laguerre-Gaussian beam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LG=lg(p,l,params)

yy=params.y;
xx=params.x;

[XX, YY]=meshgrid(xx,yy);
r=sqrt(XX.^2+YY.^2);
ZZ = sqrt(params.z^2 + r^2);

zr = (pi.*params.w0.^2)./params.lambda;
w = params.w0 * sqrt(1+(params.z/zr)^2);
R=params.z.*(1+(zr./params.z).^2);
Io = sqrt((2.*factorial(p))./(pi.*factorial(p + abs(l))));
v = 1./(w);
v1 = ((sqrt(2).*r)./w).^(abs(l));
v2 = exp(-((r.^2)./(w.^2)));
phi = atan2(YY,XX)+ pi/2;
rho = 2.*((r.^2)./(w.^2));
L=laguerre(p,abs(l),rho);
phase = exp(1i.*l.*phi + 1i * params.k * params.z);
im1=exp(1i.*params.k.*r.^2./(2.*R));
im2=exp(-1i.*(2.*p+abs(l)+1).*atan(params.z./zr));

LG=Io.*v.*v1.*v2.*L.*phase.*im1.*im2;
LG=LG/norm(LG, 'fro');

        
