load('/Volumes/GoogleDrive/Můj disk/ARTwin/InLocCIIRC_dataset/scans/B-315.ptx.mat', 'A');
%load('/Volumes/Elements/InLoc_dataset/database/scans/DUC1/DUC_scan_000.ptx.mat', 'A');
P = load_CIIRC_transformation('/Volumes/GoogleDrive/Můj disk/ARTwin/InLocCIIRC_dataset/alignments/B-315/transformations/trans_11.txt');
P = load_CIIRC_transformation('/Volumes/GoogleDrive/Můj disk/ARTwin/InLocCIIRC_dataset/alignments/B-315/transformations/trans_10.txt'); % if using shortlist GV
%P = load_WUSTL_transformation('/Volumes/Elements/InLoc_dataset/database/alignments/DUC1/transformations/DUC_trans_000.txt');
RGB = [A{5}, A{6}, A{7}]';
XYZ = [A{1}, A{2}, A{3}]';

XYZ = P * [XYZ; ones(1, length(XYZ))];
XYZ = bsxfun(@rdivide, XYZ(1:3, :), XYZ(4, :));

load('/Volumes/GoogleDrive/Můj disk/ARTwin/InLocCIIRC_dataset/outputs/PnP_dense_inlier/2.jpg/cutout_11_-60_0.pnp_dense_inlier.mat', 'P');

% taken from densePV...shortlist
P = [-0.0948407347909321,0.0231133695270234,-0.995224098971394,5.99540722308564;
    0.0205447633281221,-0.999472056158039,-0.0251698561579742,3.47452558566878;
    -0.995280436723292,-0.0228337712244152,0.0943158055023843,-2.40469821144251];

%load('/Volumes/Elements/InLoc_dataset/outputs/PnP_dense_inlier/IMG_0735.JPG/DUC_cutout_000_30_0.pnp_dense_inlier.mat', 'P');

dslevel = 8^-1;
Iq = imresize(imread('/Volumes/GoogleDrive/Můj disk/ARTwin/InLocCIIRC_dataset/query/2.jpg'), dslevel);
%Iq = imresize(imread('/Volumes/Elements/InLoc_dataset/query/iphone7/IMG_0735.JPG'), dslevel);
fl = 3172 * dslevel;
%fl = 4032*28/36 * dslevel;
K = [fl, 0, size(Iq, 2)/2.0; 0, fl, size(Iq, 1)/2.0; 0, 0, 1];

[ RGBpersp, XYZpersp ] = ht_Points2Persp( RGB, XYZ, K*P, size(Iq, 1), size(Iq, 2) );
figure(1);
imshow(RGBpersp);