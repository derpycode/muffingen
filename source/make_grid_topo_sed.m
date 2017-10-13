function [] = make_topo_sed(DMAX,KMAX,KMAX_FILTER,K1IN,K1OUT)

% set max depth
par_max_D = DMAX;
% load k1 file
seds=load(K1IN);
% calculate ocean layer depth distribution
D=-get_gdep(KMAX,0.1,par_max_D);
% extend for calculating mask
D(KMAX+2) = 1.0;
% replace shallow ocean areas with depth 'KMAX+1' (which will correspond to 0 m)
seds(find(seds(:,:)>(KMAX_FILTER-1)))=KMAX+1;
% replace land areas with depth 'KMAX+1' (which will correspond to 0 m)
seds(find(seds(:,:)>(KMAX+1)))=KMAX+1;
% truncate extended grid
seds_k1 = seds(2:end-1,2:end-1);
% save sediment K1 file for completeness
fprint_2D(seds_k1,[K1OUT '.k1.'],'%3i','%3i',true);
% replace ocean level with depth (m)
seds_D=D(seds_k1);
% save sediment depth file
fprint_2D(seds_D,[K1OUT '.dat.'],'%8.1f','%8i',true);
% replace ocean level with mask
seds_k1(find(seds_k1(:,:)<KMAX_FILTER))=KMAX+2;
% replace ocean level with depth (m)
seds_D=D(seds_k1);
% save sediment mask file
fprint_2D(seds_D,[K1OUT '.mask.'],'%3i','%3i',true);
