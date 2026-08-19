function []=generate_parameters()
% Parameters for the SEAS BP3-QD benchmark.

yr=365*24*60*60;

param.checkpointer=0;
param.output_interval=10;
param.checkpoint_interval=1000;
param.Nt=200000;
param.live_plot=false;
param.live_plot_interval=5;

% Preliminary compact-domain run. Use 100 m or finer for production
% benchmark comparisons; the benchmark target spacing is 25 m.
param.element_size=200;
param.xsize=80e3;  % +/- 40 km from the fault
param.ysize=45e3;  % 40 km frictional fault + 5 km creeping buffer
param.Nx=round(param.xsize/param.element_size)+1;
param.Ny=round(param.ysize/param.element_size)+1;

% Geometry and sense of motion. Use +1 for thrust and -1 for normal.
param.alpha=60;
param.motion_sign=-1;

% Provenance written into the SEAS output headers. The benchmark spec lists
% modeler and date as required header fields; override in run_config.json if a
% particular submission needs different attribution.
param.code_name='indnuc staggered-grid QD';
param.code_version='';
param.modeler='Meng Li';

% Homogeneous elastic medium.
param.rho=2670;
param.nu=0.25;
param.cs=3.464e3;
param.Biot=0;

% BP3 rate-and-state friction.
param.sigma0=50e6;
param.a0=0.010;
param.amax=0.025;
param.b0=0.015;
param.L0=0.008;
param.V0=1e-6;
param.f0=0.6;
param.H=15e3;
param.h=3e3;
param.Wf=40e3;

% Vp, VL and Vinit follow motion_sign and are resolved after the overlay below,
% not here -- see the note at the bottom of this file.
param.load_side_boundaries=true;
param.velocity_bracket=10;

% Time integration.
param.tfinal=1500*yr;
param.dt0=1;
param.dtmax=0.1*yr;
param.dt_growth=1.2;
param.friction_tolerance=5;

% A run folder can supply case-specific values without editing this file.
% Grid counts are recomputed after applying the overrides.
config_file=fullfile(pwd,'run_config.json');
supplied={};
if isfile(config_file)
    config=jsondecode(fileread(config_file));
    names=fieldnames(config);
    supplied=names;
    for k=1:numel(names)
        if ~strcmp(names{k},'tfinal_years')
            param.(names{k})=config.(names{k});
        end
    end
    if isfield(config,'tfinal_years')
        param.tfinal=config.tfinal_years*yr;
    end
    param.Nx=round(param.xsize/param.element_size)+1;
    param.Ny=round(param.ysize/param.element_size)+1;
end

% Driving velocities, resolved here so a run_config.json that sets motion_sign
% actually reaches them. They used to be assigned above the overlay, so a config
% asking for thrust got thrust's friction bracket and header with the default
% normal loading -- the bracket and the drive then disagree and the friction
% solve fails outright.
%
% motion_sign is the benchmark's sense of motion (+1 thrust, -1 normal), and is
% also the sign the spec assigns to slip, slip rate and shear traction. The
% solver's internal slip variable is the jump uy(+) - uy(-) set by the fault row
% in build_LH, which is the NEGATIVE of the spec's delta = uy(-) - uy(+)
% (eq 6). So the internal drive carries -motion_sign and write_bp3_outputs
% negates slip and shear stress on the way back out. Both halves are needed:
% fixing only the output sign simulates the wrong sense, fixing only the drive
% labels it wrongly. At alpha=90 the two senses are mirror images so the error
% is invisible, but for a dipping fault they are different physical problems --
% Omega+ is the hanging wall, and normal faulting needs it to move down-dip,
% uy(+) > uy(-).
internal_sign=-param.motion_sign;
for f={'Vp','VL','Vinit'}
    if ~any(strcmp(supplied,f{1}))
        param.(f{1})=internal_sign*1e-9;
    end
end

save('parameters.mat','param');
end
