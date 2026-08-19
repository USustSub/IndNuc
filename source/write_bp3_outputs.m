function write_bp3_outputs(param,xd,t,U,V,tau,sigma,theta,...
    surface_x,disp1,disp2,vel1,vel2)
% Write the ASCII products requested for the SEAS BP3-QD benchmark.
%
% Filenames follow the spec: the on-fault (fltst_dp*) and off-fault (srfst_fn*)
% time series carry NO extension, only the three slip/stress evolution profiles
% are named *.dat.  Headers carry the fields the spec does not mark optional --
% problem, code, modeler, date, element size, station location -- plus the
% per-column descriptions the server uses to interpret the file.

% Runs saved before the provenance fields existed still have to be writable.
if ~isfield(param,'code_name')
    param.code_name='indnuc staggered-grid QD';
end
if ~isfield(param,'modeler')
    param.modeler='Meng Li';
end

outdir='output_BP3_QD';
if ~exist(outdir,'dir')
    mkdir(outdir);
end

fault_stations=[0,2.5,5,7.5,10,12.5,15,17.5,20,25,30,35]*1e3;
fault_names={'000','025','050','075','100','125','150','175','200',...
    '250','300','350'};
motion_name='normal';
if param.motion_sign>0
    motion_name='thrust';
end
nt=numel(t);

for k=1:numel(fault_stations)
    [~,iy]=min(abs(xd-fault_stations(k)));
    filename=fullfile(outdir,['fltst_dp',fault_names{k}]);
    fid=fopen(filename,'w');
    write_common_header(fid,param,motion_name,nt);
    fprintf(fid,'# location=on fault, %.1f km down-dip distance\n',...
        fault_stations(k)/1e3);
    fprintf(fid,'# Column #1 = Time (s)\n');
    fprintf(fid,'# Column #2 = Slip (m)\n');
    fprintf(fid,'# Column #3 = Slip rate (log10 m/s)\n');
    fprintf(fid,'# Column #4 = Shear stress (MPa)\n');
    fprintf(fid,'# Column #5 = Normal stress (MPa)\n');
    fprintf(fid,'# Column #6 = State (log10 s)\n');
    fprintf(fid,'# The line below lists the names of the data fields\n');
    fprintf(fid,'t slip slip_rate shear_stress normal_stress state\n');
    fprintf(fid,'# Here is the time-series data.\n');
    % Slip and shear stress are negated into the spec's convention: the solver
    % carries the jump uy(+)-uy(-), the spec's delta is uy(-)-uy(+) (eq 6).
    % Slip rate and state are log10 magnitudes and normal stress is unsigned by
    % this, so only these two columns flip.
    values=[t(:),-U(iy,:)',log10(max(abs(V(iy,:)'),realmin)),...
        -tau(iy,:)'/1e6,sigma(iy,:)'/1e6,...
        log10(max(theta(iy,:)',realmin))];
    fprintf(fid,'%21.13E %14.6E %14.6E %14.6E %14.6E %14.6E\n',...
        values');
    fclose(fid);
end

surface_names={'srfst_fn-32','srfst_fn-16','srfst_fn-08','srfst_fn+00',...
    'srfst_fn-00','srfst_fn+08','srfst_fn+16','srfst_fn+32'};
% Nominal station positions from the spec; fn+00 and fn-00 are the two sides of
% the fault trace, actually sampled at the adjacent staggered nodes x=+/-dx/2.
surface_nominal=[-32e3,-16e3,-8e3,0,0,8e3,16e3,32e3];
for k=1:numel(surface_names)
    filename=fullfile(outdir,surface_names{k});
    fid=fopen(filename,'w');
    write_common_header(fid,param,motion_name,nt);
    fprintf(fid,'# location=on surface, %+g km distance off-fault\n',...
        surface_nominal(k)/1e3);
    fprintf(fid,'# sampled at x=%+.4f km\n',surface_x(k)/1e3);
    fprintf(fid,'# Column #1 = Time (s)\n');
    fprintf(fid,'# Column #2 = Displacement 1 (m)\n');
    fprintf(fid,'# Column #3 = Displacement 2 (m)\n');
    fprintf(fid,'# Column #4 = Velocity 1 (m/s)\n');
    fprintf(fid,'# Column #5 = Velocity 2 (m/s)\n');
    fprintf(fid,'# The line below lists the names of the data fields\n');
    fprintf(fid,'t disp_1 disp_2 vel_1 vel_2\n');
    fprintf(fid,'# Here is the time-series data.\n');
    values=[t(:),disp1(k,:)',disp2(k,:)',vel1(k,:)',vel2(k,:)'];
    fprintf(fid,'%21.13E %14.6E %14.6E %14.6E %14.6E\n',values');
    fclose(fid);
end

profile_mask=xd<=param.Wf;
profile_indices=find(profile_mask);
stride=max(1,round(500/param.element_size));
profile_indices=profile_indices(1:stride:end);
if profile_indices(end)~=find(profile_mask,1,'last')
    profile_indices(end+1)=find(profile_mask,1,'last');
end

% Negative scales on slip and shear stress for the same eq (6) convention flip
% applied to the on-fault time series above; normal stress is unaffected.
write_profile(fullfile(outdir,'slip.dat'),param,'slip','Slip (m)',...
    xd(profile_indices),t,V,U(profile_indices,:),-1);
write_profile(fullfile(outdir,'shear_stress.dat'),param,'shear_stress',...
    'Shear stress (MPa)',xd(profile_indices),t,V,...
    tau(profile_indices,:),-1e-6);
write_profile(fullfile(outdir,'normal_stress.dat'),param,'normal_stress',...
    'Normal stress (MPa)',xd(profile_indices),t,V,...
    sigma(profile_indices,:),1e-6);
end

function write_common_header(fid,param,motion_name,nt)
fprintf(fid,'# This is the file header:\n');
fprintf(fid,'# problem=SEAS Benchmark BP3-QD\n');
fprintf(fid,'# code=%s\n',param.code_name);
if isfield(param,'code_version') && ~isempty(param.code_version)
    fprintf(fid,'# version=%s\n',param.code_version);
end
fprintf(fid,'# modeler=%s\n',param.modeler);
fprintf(fid,'# date=%s\n',datestr(now,'yyyy/mm/dd'));
fprintf(fid,'# element size=%g m\n',param.element_size);
fprintf(fid,'# motion=%s\n',motion_name);
fprintf(fid,'# dip angle=%g degrees\n',param.alpha);
fprintf(fid,'# num time steps=%d\n',nt);
end

function write_profile(filename,param,field_name,description,xd,t,V,field,scale)
n=numel(xd);
fid=fopen(filename,'w');
fprintf(fid,'# This is the file header:\n');
fprintf(fid,'# problem=SEAS Benchmark BP3-QD\n');
fprintf(fid,'# modeler=%s\n',param.modeler);
fprintf(fid,'# date=%s\n',datestr(now,'yyyy/mm/dd'));
fprintf(fid,'# code=%s\n',param.code_name);
if isfield(param,'code_version') && ~isempty(param.code_version)
    fprintf(fid,'# code version=%s\n',param.code_version);
end
fprintf(fid,'# element size=%g m\n',param.element_size);
fprintf(fid,'# Row #1 = Distance down dip (m) with two zeros first\n');
fprintf(fid,'# Column #1 = Time (s)\n');
fprintf(fid,'# Column #2 = Max slip rate (log10 m/s)\n');
fprintf(fid,'# Columns #3-%d = %s\n',n+2,description);
fprintf(fid,['# Computational domain size: down-dip %g km, ' ...
    'distance off fault %g km, dip %g degrees\n'], ...
    param.ysize/1e3,param.xsize/2e3,param.alpha);
fprintf(fid,'# The line below lists the names of the data fields\n');
fprintf(fid,'xd\n');
fprintf(fid,'t max_slip_rate %s\n',field_name);
fprintf(fid,'# Here are the data\n');
fprintf(fid,'%14.6E ',[0,0,xd(:)']);
fprintf(fid,'\n');
for it=1:numel(t)
    fprintf(fid,'%21.13E %14.6E ',t(it),log10(max(abs(V(:,it)))));
    fprintf(fid,'%14.6E ',field(:,it)*scale);
    fprintf(fid,'\n');
end
fclose(fid);
end
