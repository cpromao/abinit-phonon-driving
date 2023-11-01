clear variables
close all

%% Inputs

filen = 'gamma_PHBST.nc';

phonon_d1 = 106;
phonon_d2 = 108;
cep1 = 0;
cep2 = 0;

qpt_gamma = 1;

t_0 = 0; %peak time in ps
t_1 = t_0 + 5.5; %final time 

phonon_c = 5;
coupling111 = 0.00199488097856232 ; % Units: eV Angstrom^-3 amu^-1.5 (M L^2 T^-2 * M^-1.5 L^-3 = M^-0.5 L^-1 T^-2) 
exportname = '5_106_108.mat';

% phonon_c = 107;
% coupling111 = 0.0561535485560740 ;
% exportname = '107_106_108.mat';

%% Constants

pm = 1822.89; % 1 amu in electron mass units
b2A = 1.88973; % 1 Angstrom in Bohr radii
a2bohr = 0.529177; %1 Angstrom in Bohr radii
T2eV = 0.00413566972; % 1 THz in eV
hbar = 1; % Hartree units
hbar_SI = 1.054571817e-34; %J/Hz
e2C = 1.60218e-19; %1 e in C
emu2k = 9.10938e-31; %1 emu in kg
amu2k = 1.66054e-27; %1 amu in kg
nmag_SI = 5.050783699e-27; %nuclear magneton in J/T
nmag = nmag_SI/1.85480201566e-23; %1 nuclear magneton in Hartree units
eV2H = 0.037; %1 eV in Ha
Heb2V = 5.14220674763e11; %1 Ha/e/bohr in V/m
H2s = 2.4188843265857e-17; % convert atomic time units to s
eV2Hz = 241799050402293; % convert eV to Hz
amu2kg = 1.66054e-27; % 1 amu in kg
J2eV = 6.242e18; %1 J in eV 
Aaps2eV = 1e-20 * amu2kg * 1e12^2 * J2eV ; % convert A^2 amu/ps^2 into eV
permittivity = 8.8541878128e-12; % e^2 eV-1 um^-1
permittivity_aaps = permittivity * Aaps2eV / 1e4;

%% Load stuff

eigend_abi = ncread(filen, 'phdispl_cart'); %Eigendisplacements
qpoints = ncread(filen, 'qpoints'); %List of qpoints
mass_amu = ncread(filen, 'atomic_mass_units'); %atomic masses in amu
species = ncread(filen, 'atom_species'); %atomic species labels
bec = ncread(filen, 'becs_cart'); %Born effective charge tensors 
phfrq = ncread(filen, 'phfreqs'); %Phonon energies in eV
rprim = ncread(filen, 'primitive_vectors');
xred = transpose(ncread(filen, 'reduced_atom_positions')); %reduced coordinates
atomic_number = ncread(filen, 'atomic_numbers');
nqpt = size(eigend_abi,4); %number of q points

natom = size(species,1); %number of atoms
nphonon = size(eigend_abi,2); %number of phonon bands

mass = mass_amu.*pm; %convert masses to emu

vol = dot(rprim(:,1),cross(rprim(:,2),rprim(:,3)));

coupling111 = coupling111/Aaps2eV; %-> A^2 amu/ps^2

%% Get cartesian coordinates

for atom = 1:natom
    xcart(atom,:) = xred(atom,1)*rprim(1,:)+xred(atom,2)*rprim(2,:)+xred(atom,3)*rprim(3,:);
end

%% Reformat eigendisplacements (bohr) and get unit eigenvectors (unitless)

for phonon = 1:nphonon
	for atom = 1:natom
		amf(atom) = sqrt((mass(species(atom))));
        for qpt = 1:nqpt
            for dir = 1:3
                eigendr(atom,dir,phonon,qpt) = eigend_abi(1,(atom-1)*3+dir,phonon,qpt)/a2bohr; %real part
                eigendi(atom,dir,phonon,qpt) = eigend_abi(2,(atom-1)*3+dir,phonon,qpt)/a2bohr; %imaginary part
                eigend(atom,dir,phonon,qpt) = complex(eigendr(atom,dir,phonon),eigendi(atom,dir,phonon));
                eivecr(atom,dir,phonon,qpt) = eigendr(atom,dir,phonon)*amf(atom);
                eiveci(atom,dir,phonon,qpt) = eigendi(atom,dir,phonon)*amf(atom);
                eivec(atom,dir,phonon,qpt) = complex(eivecr(atom,dir,phonon,qpt),eiveci(atom,dir,phonon,qpt));
            end
            eivec_n(atom,:,phonon,qpt) = eivec(atom,:,phonon,qpt)/norm(eivec(atom,:,phonon,qpt));
            eigend_n(atom,:,phonon,qpt) = eigend(atom,:,phonon,qpt)/norm(eivec(atom,:,phonon,qpt));
        end
    end
end

%% Calculate mode effective charges (e / √emu)

mecc = zeros(3,nphonon,nqpt);
mecc2 = zeros(3,nphonon,nqpt);

for qpt = 1:nqpt
    for phonon = 1:nphonon
        for atom = 1:natom
            amf(atom) = sqrt((mass(species(atom))));
            mecc(:,phonon,qpt) = mecc(:,phonon,qpt) + bec(:,:,atom)*transpose(eigend(atom,:,phonon,qpt))/amf(atom);
            mecc2(:,phonon,qpt) = mecc2(:,phonon,qpt) + bec(:,:,atom)*transpose(eigend(atom,:,phonon,qpt));
            eigend_norm(atom,phonon,qpt) = dot(transpose(eigend(atom,:,phonon,qpt)),eigend(atom,:,phonon,qpt));
        end
        mecc_e(:,phonon,qpt) = mecc(:,phonon,qpt)/sqrt(sum(eigend_norm(:,phonon,qpt),1));
        mecc_e2(:,phonon,qpt) = mecc2(:,phonon,qpt)/sqrt(sum(eigend_norm(:,phonon,qpt),1));
        mec_scalar(phonon,qpt) = norm(mecc_e(:,phonon,qpt));
        mec_scalar2(phonon,qpt) = norm(mecc_e2(:,phonon,qpt));
        mecc_n(:,phonon,qpt) = mecc_e(:,phonon,qpt)/norm(mecc_e(:,phonon,qpt));
    end
end

%% Solve eq 3

% Q units: Angst √amu (L M^0.5)
% Z Q E = energy 
% Z L M^0.5 C^-1 L^-1 = 1
% Z = C M^-0.5


omega_d1 = phfrq(phonon_d1,qpt_gamma)*eV2Hz /1e12;
omega_d2 = phfrq(phonon_d2,qpt_gamma)*eV2Hz /1e12;
omega_c = phfrq(phonon_c,qpt_gamma)*eV2Hz /1e12;

drivefreq = mean([omega_d1,omega_d2]);

linew_d1 = 0.1*omega_d1;
linew_d2 = 0.1*omega_d2;
linew_c = 0.1*omega_c;

mec_d1 = dot(mecc_e(:,phonon_d1,qpt_gamma), [0, 1, 0]); % mec in e / √emu
mec_d_SI1 = mec_d1*e2C/emu2k^0.5; % mec in C / √kg

mec_d2 = dot(mecc_e(:,phonon_d2,qpt_gamma), [0, 0, 1]); 
mec_d_SI2 = mec_d2*e2C/emu2k^0.5; 

E_0 = 30 / sqrt(2) * 1e6 / 1e-2; % V/m = eV/e/m
E_au = E_0/Heb2V; %Ha/e/bohr

prefactor_SI1 = mec_d_SI1*E_0; % V/m C √kg = J/C/m C /√kg = J/m /√kg = m √kg Hz^2 (M^0.5 L T^-2)
prefactor1 = prefactor_SI1*amu2kg^-0.5*1e10 / 1e12^2; % Angstrom √amu THz^2

prefactor_SI2 = mec_d_SI2*E_0; 
prefactor2 = prefactor_SI2*amu2kg^-0.5*1e10 / 1e12^2; 

tau = 0.1;
const = sqrt(8*log(2));
time = linspace(0,t_1,1000);


syms Q_d1(t) Q_d2(t) Q_c(t) 

ode1 = diff(Q_d1(t),2) ==  prefactor1*exp((-(t-t_0)^2)/(2*(tau/const)^2))*cos(drivefreq*t*2*pi+cep1) - linew_d1*diff(Q_d1(t)) - (omega_d1*2*pi)^2*Q_d1(t) -coupling111*Q_d2(t)*Q_c(t);
ode2 = diff(Q_d2(t),2) ==  prefactor2*exp((-(t-t_0)^2)/(2*(tau/const)^2))*cos(drivefreq*t*2*pi+cep2) - linew_d2*diff(Q_d2(t)) - (omega_d2*2*pi)^2*Q_d2(t) -coupling111*Q_d1(t)*Q_c(t); 
ode3 = diff(Q_c(t),2) == -linew_c*diff(Q_c(t)) -(omega_c*2*pi)^2*Q_c(t) -coupling111*Q_d2(t)*Q_d1(t);

odes = [ode1; ode2; ode3];
vars = [Q_d1(t), Q_d2(t) Q_c(t)];
[newEqs, newVars] = reduceDifferentialOrder(odes, vars);
[M,F] = massMatrixForm(newEqs, newVars);
M = odeFunction(M, newVars);
F = odeFunction(F, newVars);
opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'mass',M,'OutputFcn',@odeplot); %ODE solver settings
sol = ode15s(F,time,[0 0 0 0 0 0], opts);

%% Find smallest interatomic distance for each atom

distances = zeros(natom,natom,3,3,3);

for atom1 = 1:natom
    for atom2 = 1:natom
        for i = 1:3
            for j = 1:3
                for k = 1:3
                    distances(atom1,atom2,i,j,k) = norm(xcart(atom1,:)-xcart(atom2,:)+rprim(:,1)*(i ~= 1)*(-2*(i ~= 2)+1)+rprim(:,2)*(j ~= 1)*(-2*(j ~= 2)+1)+rprim(:,3)*(k ~= 1)*(-2*(k ~= 2)+1));
                    if distances(atom1,atom2,i,j,k) == 0
                        distances(atom1,atom2,i,j,k) = NaN;
                    end
                end
            end
        end
    end
end

mindistances = min(distances,[],[3 4 5]);
mindistance = min(distances,[],[1 2 3 4 5]);

%% Evaluate Lindemman criterion 

Qtest = max([max(sol.y(1,:)),max(sol.y(2,:))]);
Qtest_au = Qtest / b2A *pm^0.5;

for atom = 1:natom
    for dir = 1:3
        disp{atom,dir} = Qtest*eivec(atom,dir,phonon_d1,qpt_gamma)./(amf(atom));
    end
    displ{atom} = norm([disp{atom,1},disp{atom,2},disp{atom,3}]);
    rmsdispatom(atom) = rms((displ{atom}));
end

for atom = 1:natom
    lindemann_atom(atom) = rmsdispatom(atom)/min(mindistances(atom,:))*100;
end

lindemann = max(lindemann_atom(:)); %Lindemann factor

% S = movmean(sol.y(3,:),[5000,500]);
S = movmean(sol.y(3,:),[1000,500]);

%% Solve Eq 10

freqlim = 120;
freq1 = linspace(0,freqlim,100000);

for i = 1:size(freq1,2)
    freq2 = omega_d1;
    chi_c(i) = - (coupling111) * (mec_d1 *mec_d2 * pm) / ((omega_c^2 - freq1(i)^2 + 1i*freq1(i)*linew_c/(2*pi)) * (omega_d1^2 - (freq1(i)-freq2)^2 + 1i*(freq1(i)-freq2)*linew_d1/(2*pi)) * (omega_d2^2 - freq2^2 + 1i*(freq2)*linew_d2/(2*pi))*(2*pi)^6);
end

chi_c_SI = chi_c * (e2C^2 * 1e-24 * amu2k^-2 * 1e20);% *(1e10 * amu2k^0.5); % convert to SI units

%% FFT

t_eval = linspace(0,5,5000);
Q_eval = deval(sol,t_eval);

fft_Qc = fft(Q_eval(3,:),1000);

%% Save variables for external plotting

i = phonon_c;

times{i} = sol.x;
Q_ci{i} = sol.y(3,:);
chi_ci{i} = chi_c;
freqi{i} =  freq1;
Si{i} = S;
ffti{i} = fft_Qc;

save(exportname,"times","Q_ci","chi_ci","freqi","Si","ffti")

%% Make figures

linew = 1;
time_append = linspace(-2,0,10);
Q_cappend = linspace(0,0,10);
plottime = cat(2,time_append,sol.x)*1000;
Q_cplot = cat(2,Q_cappend,sol.y(3,:));
Splot = cat(2,Q_cappend,S);
orange = '#D95319';
lblue = '#1597EE';
lime = '#32CD32';
scale = 1e3;

close all 

% figure
% plot(sol.x,sol.y(1,:))

% figure
% plot(sol.x,sol.y(2,:))

f1 = figure;

hold on

a = area(plottime,Splot*scale,'EdgeColor',"none",'FaceColor',orange);
a.FaceAlpha = 0.5;

a2 = area(plottime,-Splot*scale,'EdgeColor',"none",'FaceColor',lblue);
a2.FaceAlpha = 0.5;

p = plot(plottime,Q_cplot*scale,'LineWidth',linew,'Color',orange);
p2 = plot(plottime,-Q_cplot*scale,'LineWidth',linew,'Color',lblue);

xlabel('t [fs]')
ylabel('Q_c [Å amu^{1/2} \times 10^{-3}]')
xlim([-300,1100])

% ylim([-0.2,0.2])

ax = gca;

box on

hold off

f2 = figure;

tiledlayout(2,1,"TileSpacing","compact")

nexttile

plot(abs(fft_Qc/max(abs(fft_Qc))),'LineWidth',linew,'Color','blue');

ylabel('||FFT(\omega)||')
ax1 = gca;
set(ax1,'XTickLabel',[])

xlim([0,110])

nexttile

plot(freq1,log10(abs(chi_c_SI)),'LineWidth',linew,'Color',lime);

ylabel('|\chi| [(m/V)^2]')
xlabel('\omega / 2\pi [THz]')
xlim([0,110])


yticklabels({'10^{0}','10^{2}','10^{4}'})

if phonon_c == 5
    ylim([-0.5,5.5])
else
    ylim([0,4])
end


ax2 = gca;


%% Do some post-adjustment of figures

xt2 = [0 25 50 75 100];
xt1 = xt2;

set(ax1,'xtick',xt1)
set(ax2,'xtick',xt2)

xt = [-250 0 250 500 750 1000];
set(ax,'xtick',xt)

if phonon_c == 5
     yt = [0 2 4];
     set(ax2,'ytick',yt)
    yt2 = [-1.5 -1 -0.5 0 0.5 1 1.5];
    set(ax,'ytick',yt2)
else
    yt = [0 2 4];
    set(ax2,'ytick',yt)
    yt2 = [-0.2 -0.1 0 0.1 0.2];
    set(ax,'ytick',yt2)
end

ax2.XAxis.FontSize = 12;
ax2.XLabel.FontSize = 16;
ax2.YAxis.FontSize = 12;
ax2.YLabel.FontSize = 16;
ax1.YAxis.FontSize = 12;
ax1.YLabel.FontSize = 16;
ax.XAxis.FontSize = 12;
ax.XLabel.FontSize = 16;
ax.YAxis.FontSize = 12;
ax.YLabel.FontSize = 16;
