%% Authors: H. Nguyen-Xuan, Khanh Nguyen Chau, Khai Nguyen Chau
%  Email:ngx.hung@hutech.edu.vn
%  Codes derived from [Polytopal composite finite elements, Computer
%  Methods in Applied Mechanics and Engineering, 355, 405-437, 2019]

clear all
close all
clc

format long

tic

addpath(genpath('cantilever3DSundar'))
addpath(genpath('functions'))

% load cant3DmeshAP10elem.mat
% load cant3DmeshAP50elem.mat
load cant3DmeshAP100elem.mat

Face = Elements.face.vertex_indices;
NELE = size(Face, 1);

for el = 1:NELE
    Face{el} = Face{el} + 1;
end

for el = 1:numel(Elements.face.vertex_indices)
    Elements.face.vertex_indices{el} = Elements.face.vertex_indices{el} + 1;
end

for fa = 1 : numel(Elements.cell.face_indices)
    Elements.cell.face_indices{fa} = Elements.cell.face_indices{fa} + 1;
end
node(:, 1) = node(:, 1) - 1;
node(:, 2) = node(:, 2) - 1;

Mesh.vertices = node;
Mesh.nvertex = size(node, 1);
Mesh.faces = Elements.face.vertex_indices;
Mesh.nface = numel(Elements.face.vertex_indices);
% for i = 1 : numel(element)
%     Mesh.ele{i} = {element{i}, Elements.cell.face_indices{i}};
% end
% Mesh.bdryfidx = [];
% Mesh.nele = numel(element);
% Mesh.eleidx = [];
% Mesh.nBf = [];

numnode = size(node, 1);
numelem = numel(element);
ndof = 3;

cenDOF = 'no';
flMode = 'PCEn';
% flMode = 'PEn';

% trifg = 'fan'; % specify triangulation method
trifg = 'mid';

% flags to print node numbers and element numbers
pflag.pMesh = 'yes';
pflag.pNode = 'no';
pflag.pElem = 'no';

if strcmp(pflag.pMesh, 'yes')
    figure(1)
    Element = Face(1:NELE)'; %Only plot the first block
    MaxNVer = max(cellfun(@numel, Element)); %Max. num. of vertices in mesh
    PadWNaN = @(E) [E NaN(1, MaxNVer - numel(E))]; %Pad cells with NaN
    ElemMat = cellfun(PadWNaN, Element, 'UniformOutput', false);
    ElemMat = vertcat(ElemMat{:}); %Create padded element matrix
    aux_node = node;
    tmp_node = node(:, 3);
    aux_node(:, 3) = node(:, 2);
    aux_node(:, 2) = tmp_node;
    patch('Faces', ElemMat, 'Vertices', aux_node, 'FaceVertexCData', hsv(1), 'FaceColor', 'r'); %pause(1e-6)
    alpha(0.1); view(3); axis equal
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    viewAz = 45;
    vieEl = 25;
    view(viewAz, vieEl);
end

if strcmp(pflag.pNode, 'yes')
    % plot node numbers
    for in = 1:numnode
        xc = node(in, 1) + 0.001;
        yc = node(in, 2) - 0.001;
        zc = node(in, 3) - 0.001;
        text(xc, yc, zc, num2str(in), 'color', 'blue');
    end
end

if strcmp(pflag.pElem, 'yes')
    %print element numbers
    for iel = 1:numelem
        econ = element{iel};
        nod = node(econ, :);
        text(mean(nod(:, 1)), mean(nod(:, 2)), mean(nod(:, 3)), num2str(iel), 'color', 'red');
    end
end

L = 5;
a = 1;
b = 1;
clear bnodes
% get boundary nodes
bnodes = find(node(:, 3) > L - 1e-3 & node(:, 3) < L + 1e-3);

% % % % loop over the Face...
bface={};
for in = 1:size(Face,1)
    ele_con=Face{in};
    cood=node(ele_con,:);
    trace_face=find(cood(:,3)==0);
    if length(trace_face)>=3
        %disp('trace_face')
        %elebot=[elebot; fliplr(ele_con)];
        bface=[bface; (ele_con)];
    end
end

% material properties...
E = 3e3;
% nu = 0.3;
nu = 0.4999999;

basisLAB = 'PL';
% basisLAB = 'Wachspress';
% basisLAB = 'CorrectedWachspress';

% total system degrees of freedom
totalUnknown = ndof * numnode;

num_shp = cellfun(@numel, element);
num_shp_max = max(num_shp);

el_conn = zeros(num_shp_max, numelem);
for i = 1 : numelem
    el_conn(1 : num_shp(i), i) = element{i};
end

% initialize the system matrices...
F = zeros(totalUnknown, 1);

force = 1;

[P, ~, W, ~] = getQuadData(3);
% loop over boundary elements
for bf = 1:numel(bface)
    cface = bface{bf};
    cface_coord = node(cface, :);
    Coord = [cface_coord; mean(cface_coord)];
    sctr = cface;
    nnel = length(sctr);
    
    sctry = ndof .* sctr - 1;    
    [tri, tri_coord] = getSurfaceTriangulation(cface, cface_coord, trifg);
    
    for it = 1:size(tri, 1)
        ctri = tri(it, :);
        ctri_coord = tri_coord(ctri, :);
        
        for igp = 1 : numel(W)
            [NT3, dNdsT3] = T3ShapeFnc(P(igp, :));
            xy = NT3*ctri_coord;
            J = dNdsT3 * ctri_coord;
            dNdxT3 = J \ dNdsT3;
            
            % hold on
            % plot3(xy(1), xy(2), xy(3), '*')
            
            switch trifg
                case 'fan'
                    NPL = zeros(nnel, 1);
                    NPL(ctri) = NT3;
                case 'mid'
                    NPL = zeros(nnel, 1);
                    NPL(ctri(2:3)) = NT3(2:3);
                    NPL = NPL + repmat(NT3(1) / nnel, nnel, 1);
                otherwise
                    error('Not implemented yet!')
            end
            
            % dNdx = zeros(nnel, 3);
            % dNdx(ctri(2:3), :) = dNdxT3(:, 2:3)';
            % dNdx = dNdx + repmat(dNdxT3(:, 1)' / nnel, nnel, 1);
            
            jac_det = norm(cross(J(1, :), J(2, :)));
            
            F(sctry) = F(sctry) - NPL * force/(2*a*2*b) * jac_det * W(igp);
        end
    end
end

K = computePCEnStiffnessMatrix3D(E, nu, element, node, Elements, basisLAB, flMode, trifg, cenDOF);

% get the boundary conditions....
bcdof = [ndof * bnodes - 2; ndof * bnodes - 1; ndof * bnodes];
xpt = node(bnodes, 1);
ypt = node(bnodes, 2);
zpt = node(bnodes, 3);
I = 4*a*b^3 / 3;
solu = @(x, y, z) PrismaticBeamExactSolu(a, b, E, I, nu, force, x, y, z);
num_bnodes = numel(bnodes);
UX = zeros(num_bnodes, 1);
UY = zeros(num_bnodes, 1);
UZ = zeros(num_bnodes, 1);

for i = 1 : num_bnodes
    % hold on
    %     plot3(xpt(i), ypt(i), zpt(i), '*')
    aux_solu = solu(xpt(i), ypt(i), zpt(i));
    UX(i) = aux_solu.ux;
    UY(i) = aux_solu.uy;
    UZ(i) = aux_solu.uz;
end
bcval = [UX; UY; UZ];

% SOLVE SYSTEM
disp([num2str(toc), '   SOLVING SYSTEM'])
freedof = setdiff((1:ndof * numnode)', bcdof);
U = zeros(ndof * numnode, 1);
U(bcdof) = bcval;
F(freedof) = F(freedof) - K(freedof, bcdof) * bcval;

% solving system of equations for free dofs
U(freedof) = K(freedof, freedof) \ F(freedof);

disp(['Num DOFs = ', num2str(ndof * numnode - numel(bcdof))])

StrainEnergy = 0.5*F'*U;

uxh = U(1:ndof:ndof * numnode);
uyh = U(2:ndof:ndof * numnode);
uzh = U(3:ndof:ndof * numnode);

ux_ana = zeros(numnode, 1);
uy_ana = zeros(numnode, 1);
uz_ana = zeros(numnode, 1);
for i = 1 : numnode
    aux_solu = solu(node(i, 1), node(i, 2), node(i, 3));
    ux_ana(i) = aux_solu.ux;
    uy_ana(i) = aux_solu.uy;
    uz_ana(i) = aux_solu.uz;
end
% U(1:ndof:ndof * numnode) = uxhh;
% U(2:ndof:ndof * numnode) = uyhh;
% U(3:ndof:ndof * numnode) = uzhh;

% compute the L2 error...
%%%%%%% at nodes
Err = 0; De = 0;

for i = 1:numnode
    Err = Err + ((ux_ana(i, 1) - uxh(i))^2 + (uy_ana(i, 1) - uyh(i))^2 + (uz_ana(i, 1) - uzh(i))^2);
    De = De + (ux_ana(i, 1)^2 + uy_ana(i, 1)^2 + uz_ana(i, 1)^2);
end

Relerrdisp = sqrt(Err / De)
%==========================================================================
%               END PROCESSING/SOLUTION PHASE
%==========================================================================

[L2, H1] = computeErrorNormPCEn3D(E, nu, element, node, Elements, U, solu, basisLAB, flMode, trifg, 'no');
disp(['Displacement Norm = ', num2str(L2, 20)]);
disp(['Energy Norm = ', num2str(H1, 20)]);

% [L2, H1] = computeErrorNormLS43D(element, node, Elements, ndof, U, solu, matmtx);
% disp(['Displacement Norm = ', num2str(L2, 20)]);
% disp(['Energy Norm = ', num2str(H1, 20)]);

exportPCEn3DDisplToVTU(node, element, Elements.face.vertex_indices, Elements.cell.face_indices, U, 'PrismaticBeamDisplSundarMesh')