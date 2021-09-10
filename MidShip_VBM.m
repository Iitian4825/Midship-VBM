function body = read_input_file(body,inp_fname)
body.npan = npan;
body.pan = pan;
function body = calc_grid_properties(body)

% Required inputs

npan = body.npan;

% Your code here to compute centroids, normals and areas of each panel
%----------------------------------------------------------------
for i = 1:npan;
    panel=body.pan(:,:,i);
 
    d1=norm(panel(:,2)-panel(:,1));
    d2=norm(panel(:,3)-panel(:,2));
    d3=norm(panel(:,4)-panel(:,3));
    d4=norm(panel(:,1)-panel(:,4));
 
    cen(i,:)=(d1*(panel(:,1)+panel(:,2))+d2*(panel(:,2)+panel(:,3))+d3*(panel(:,3)+panel(:,4))+d4*(panel(:,4)+panel(:,1)))/(2*(d1+d2+d3+d4));
 
    D1=(panel(:,3)-panel(:,1)).';
    D2=(panel(:,4)-panel(:,2)).';
 
    normal=cross(D1,D2);
    n(i,:)=normal./norm(normal);
 
    area(i,1)=norm(normal)/2;
end
%----------------------------------------------------------------

% Added additional output fields to output structure

body.cen = cen;
body.nor = n;
body.area = area;
function body = hydrostatics_check(body)

% Required inputs

xbody = body.xbody;
vcg = body.vcg;
area = body.area;
n = body.nor;
cen = body.cen;

% Your code here to compute centroids, normals and areas of each panel
%----------------------------------------------------------------
Awp = sum(area.*n(:,3));
dispx = -sum(area.*n(:,1).*cen(:,1));
dispy = -sum(area.*n(:,2).*cen(:,2));
dispz = -sum(area.*n(:,3).*(cen(:,3)+xbody(3)));
xb = -((sum(area.*n(:,3).*(cen(:,3)+xbody(3)).*cen(:,1)))./dispz);
yb = -((sum(area.*n(:,3).*(cen(:,3)+xbody(3)).*cen(:,2)))./dispz);
zb = -((sum(area.*n(:,3).*(cen(:,3)+xbody(3)).*(cen(:,3)/2+xbody(3)/2-xbody(3))))./dispz);
cob=[xb yb zb];
Ix = sum(area.*n(:,3).*cen(:,2).*cen(:,2));
BM = Ix./dispz;
GM = BM+zb-vcg;
%----------------------------------------------------------------

% Added additional output fields to output structure

body.Awp = Awp;
body.cob = cob;
body.dispx = dispx;
body.dispy = dispy;
body.dispz = dispz;
body.GM = GM;
function body = weight_distribution(body)

rho = body.rho;
g = body.g;
Lpp = body.Lpp;
disp = body.dispz;

a = 0.2*Lpp;
b = 0.3045*Lpp;

% Longitudinal extent of the vessel

xmin = min(min(body.pan(1,:,:)));
xmax = max(max(body.pan(1,:,:)));

xdist = linspace(xmin,xmax,20);

% Export outputs

body.xdist = xdist;
body.wdist = wdist;
function VBM = calculate_vbm_still(body)

% Constants

rho = body.rho;
g = body.g;

% Required Inputs

xbody = body.xbody;
xdist = body.xdist;
wdist = body.wdist;

pan_cen = body.cen;
pan_area = body.area;
pan_nz = body.nor(:,3);

% Global coordinates of the panel centroids

X = pan_cen(:,1) + xbody(1);
Z = pan_cen(:,3) + xbody(3);

% Longitudinal spacing at which VBM is evaluated

dxdist = mean(diff(xdist));

for i = 1:20
   xk = [xdist(i) - dxdist/2,xdist(i) + dxdist/2];
   x = zeros(1,length(X));
   z = zeros(1,length(X));
   for k = 1:length(X)
       if pan_cen(k,1)>xk(1) && pan_cen(k,1)<xk(2)
           x(k) = X(k,1);
       end
   end
   for j = 1:length(x)
       if x(j) ~= 0
           z(j) = Z(j,1);
       end
   end
   s = sum(-rho*g*z'.*pan_area.*pan_nz);
   if i == 1 || i == 20
       bdist(i) = (2/dxdist)*s;
   else
       bdist(i) = s/dxdist;
   end
end
VBM =cumtrapz(xdist, cumtrapz(xdist,(wdist-bdist)))
%--------------------------------------------------------------------
function body = calculate_RAO(body)
 
% Constants
 
rho = body.rho;
g = body.g;
 
% Required inputs
 
xbody = body.xbody;
xprdct = body.xprdct; kyy = xprdct(2,2);    % Pitch radius of gyration
w = body.w;
bet = body.dir * pi / 180;
 
zeta = body.zeta;
 
pan_cen = body.cen;
pan_area = body.area;
pan_nx = body.nor(:,1);
pan_nz = body.nor(:,3);
 
xB = body.cob(1);
xG = xB;
 
% Global coordinates of the panel centroids
 
X = pan_cen(:,1) + xbody(1);
Y = pan_cen(:,2) + xbody(2);
Z = pan_cen(:,3) + xbody(3);
 
% Wave number
 
k = w.^2/g;
 
F3 = zeros(1,numel(w));
F5 = zeros(1,numel(w));
xi = zeros(2,numel(w));
 
 
% Your code here--------------------------------------------------------------------
 
M = zeros(2,2);
C = zeros(2,2);
B = zeros(2,2);
 
% Mass matrix:
 
M(1,1) = -rho*sum(pan_area(:,1).*pan_nz(:,1).*Z(:,1))
M(1,2) = rho*sum(pan_area(:,1).*pan_nz(:,1).*Z(:,1).*pan_cen(:,1))
M(2,1) = rho*sum(pan_area(:,1).*pan_nz(:,1).*Z(:,1).*pan_cen(:,1))
M(2,2) = -rho*sum(pan_area(:,1).*pan_nz(:,1).*Z(:,1).*(kyy^2))
 
% Stiffness matrix:
 
C(1,1) = rho*g*sum(pan_area(:,1).*pan_nz(:,1))
C(1,2) = -rho*g*sum(pan_area(:,1).*pan_nz(:,1).*pan_cen(:,1))
C(2,1) = -rho*g*sum(pan_area(:,1).*pan_nz(:,1).*pan_cen(:,1))
C(2,2) = rho*g*sum(pan_area(:,1).*pan_nz(:,1).*(pan_cen(:,1).^2))
 
% Damping matrix:
 
B(1,1) = 0.6*sqrt(M(1,1)*C(1,1))
B(2,2) = 0.6*sqrt(M(2,2)*C(2,2)) 
for j = 1:numel(k)
    F3(1,j) = rho*g*sum(pan_area.*pan_nz.*exp(k(j).*Z(:,1)).*exp(-1i*(k(j).*X(:,1).*cos(bet)+k(j).*Y(:,1).*sin(bet))))
    F5(1,j) = rho*g*sum(pan_area.*((pan_cen(:,3).*pan_nx(:,1))-(pan_cen(:,1).*pan_nz(:,1))).*exp(k(j).*Z).*exp(-(1i.*k(j).*X.*cos(bet)+1i.*k(j).*Y.*sin(bet))))
end
F = [F3;F5];
eta3a = zeros(1,numel(w));
eta5a = zeros(1,numel(w));
eta = [eta3a;eta5a];
for j = 1:numel(w)
    eta(:,j) = inv(-(M*w(j)^(2))+(B*1i*w(j))+C)*F(:,j)
end
xi(1,:)=eta(1,:);
xi(2,:)=eta(2,:);
 
 
 
%--------------------------------------------------------------------
 
% Export outputs
 
body.xi3 = xi(1,:);
body.xi5 = xi(2,:);
function VBM = calculate_vbm_wave(body)
 
% Constants
 
rho = body.rho;
g = body.g;
 
% Required Inputs
 
xbody = body.xbody;
xprdct = body.xprdct;
xdist = body.xdist;
zeta = body.zeta;
 
pan_cen = body.cen;
pan_area = body.area;
pan_nx = body.nor(:,1);
pan_nz = body.nor(:,3);
 
w = body.w;
xi3 = body.xi3;
xi5 = body.xi5;
bet = body.dir * pi / 180;
 
% Global coordinates of the panel centroids
 
X = pan_cen(:,1) + xbody(1);
Y = pan_cen(:,2) + xbody(2);
Z = pan_cen(:,3) + xbody(3);
 
% Waterplane area
 
awp = sum(pan_area .* pan_nz);
 
% Wave number
 
k = w.^2/g;
 
% Longitudinal spacing at which VBM is evaluated
 
dxdist = mean(diff(xdist));
 
% Pitch radius of gyration
 
kyy = xprdct(2,2);
 
VBM = zeros(numel(w),numel(xdist));
 
% Your code here--------------------------------------------------------------------

M(2,2) = -rho*sum(pan_area(:,1).*pan_nz(:,1).*Z(:,1).*(kyy^2))
C(2,2) = rho*g*sum( pan_area(:,1).*pan_nz(:,1).*(pan_cen(:,1).^2))
B(2,2) = 0.6*sqrt(M(2,2)*C(2,2))
for i = 1:20    
    ind = find(X<(xdist(i)+dxdist/2) & X>=(xdist(i)-dxdist/2))
    for j = 1:numel(w)
        v1 = rho*g*(pan_area(ind,1).*((pan_cen(ind,3).*pan_nx(ind,1))-(pan_cen(ind,1).*pan_nz(ind,1))).*exp(k(j).*Z(ind)).*exp(-(1i.*k(j).*X(ind).*cos(bet)+1i.*k(j).*Y(ind).*sin(bet))))
        v2 = (pan_area(ind,1).*pan_nz(ind,1)).*(rho.*g.*(pan_cen(ind,1).^(2))-w(j)^(2).*(-rho.*Z(ind)*kyy^(2))+1i*w(j)*B(2,2)/awp).*xi5(j)
        v3 = (pan_area(ind,1).*pan_nz(ind,1)).*(-rho*g.*pan_cen(ind,1)).*xi3(j)
        m(j,i)=(1/dxdist)*sum(v1-v2-v3) ;    
    if i==1 | i==20
        m(j,i) = 2*m(j,i)
    end
    end
end
VBM = -cumtrapz(xdist,m,2)
end
%-------------------------
function S = wave_spectrum(w,Hs,Tp)

g = 9.80665;
function Sr = response_spectrum(S,vbm_tf)
function vbm = response_time_history(tim,wav,w,vbm_tf)

WAV = fft(wav);

[a,b] = size(tim);
Tmax = tim(end);

k = 1:1:b;
wk = (k - ones(1,b)).*((2*pi)/Tmax);
vbm_tf_new = zeros(1,b);
for i = 1:b
    if wk(i) < w(1) | wk(i) > w(end)
        vbm_tf_new(i) = 0;
    elseif wk(i) == w(1)
        vbm_tf_new(i) = vbm_tf(1);
    elseif wk(i) == w(end)
        vbm_tf_new(i) = vbm_tf(end);
    else
        ind1 = find(w < wk(i));
        ind1 = ind1(end);
        ind2 = find(w > wk(i));
        ind2 = ind2(1);
        vbm_tf_new(i) = (wk(i) - w(ind1))*(vbm_tf(ind2) - vbm_tf(ind1))/(w(ind2) - w(ind1)) + vbm_tf(ind1);
    end
end

VBM = WAV.*vbm_tf_new;
j = 1;
for i = (1 + b)/2+1:b
    VBM(i) = conj(VBM((end + 1)/2 - j));    
    j = j + 1;
end

vbm = ifft(VBM,'symmetric');


