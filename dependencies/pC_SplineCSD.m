function [new_CSD_matrix, interp_el_pos] = pC_SplineCSD(LFP, el_pos)

%% some variables for CSD
gauss_sigma = 0.00008;
cond = 0.3;
cond_top = 0.3;
diam = 0.0005;
filter_range = 5*gauss_sigma; % numeric filter must be finite in extent
el_pos = (el_pos  * 1e-6) + 1e-6; %should not be 0 at first value. Rescale depth values because code wants meters for some reason.

% compute spline iCSD:
Fcs = pC_F_cubic_spline(el_pos,diam,cond,cond_top,1e-6);
[zs,CSD_cs] = pC_make_cubic_splines(el_pos,LFP,Fcs);
[interp_el_pos,CSD_cs]=gaussian_filtering(zs,CSD_cs,gauss_sigma,filter_range);

%% matrix transformation
unit_scale = 1e-3; % A/m^3 -> muA/mm^3
CSD_matrix = CSD_cs*unit_scale;

npoints = 200; % number of points to plot in the vertical direction
le = length(interp_el_pos);
first_z = interp_el_pos(1)-(interp_el_pos(2)-interp_el_pos(1))/2; %plot starts at z1-h/2;
last_z = interp_el_pos(le)+(interp_el_pos(le)-interp_el_pos(le-1))/2; %ends at zN+h/2;
zs = first_z:(last_z-first_z)/npoints:last_z;
interp_el_pos(le+1) = interp_el_pos(le)+(interp_el_pos(le)-interp_el_pos(le-1)); % need this in for loop
j=1; %counter
for i=1:length(zs) % all new positions
    if zs(i)>(interp_el_pos(j)+(interp_el_pos(j+1)-interp_el_pos(j))/2) % > interp_el_pos(j) + h/2
        j = min(j+1,le);
    end
    new_CSD_matrix(i,:)=CSD_matrix(j,:);
end
