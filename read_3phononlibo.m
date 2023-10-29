%Requires polyfitn toolbox

clear variables

gridsize = 9; % Grid size
scale = 2; %Max scale factor in Angstrom  sqrt(amu)
phonons = [107, 106, 108];
qpt = 1;
filen = 'gamma_PHBST.nc';


amu2kg = 1.66054e-27; % 1 amu in kg
J2eV = 6.242e18; %1 J in eV 
eV2Hz = 241799050402293; % convert eV to Hz

phfrq = ncread(filen, 'phfreqs'); %Phonon energies in eV
phfrq = phfrq.*eV2Hz;
convert = 1e-20 * amu2kg * J2eV; % convert A^2 amu/s into eV

Q_line = linspace(-2+0.000001,scale,gridsize+0.000001);

load E.mat % frozen phonon energies in Ha
load modelterms.mat

E = E_107;
E = (E - E((gridsize^3+1)/2))*27.211369917461; %convert to eV

dataset = 1;
[E_grid, dataset_index] = deal(zeros(gridsize,gridsize,gridsize));
Q = zeros(gridsize,3);
E_mod = zeros(size(E));

for i = 1:gridsize
    for j = 1:gridsize
        for k = 1: gridsize
            dataset_index(i,j,k) = dataset;
            Q(dataset,:) = Q_line(transpose([i,j,k]));
            E_mod(dataset) = convert*(phfrq(phonons(1),qpt)^2/2*Q(dataset,1)^2 + phfrq(phonons(2),qpt)^2/2*Q(dataset,2)^2 + phfrq(phonons(3),qpt)^2/2*Q(dataset,3)^2);
            dataset = dataset + 1;
        end
    end
end

E = E - E_mod;

E = E(1:gridsize^3);
Q = Q(1:gridsize^3,:);


p = polyfitn(Q,E,modelterms6);

coupling = p.Coefficients(48);

% 

