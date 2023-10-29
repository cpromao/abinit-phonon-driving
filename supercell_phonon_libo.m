clear variables

%Options to adjust if necessary

filen = 'gamma_PHBST.nc';
phonon_index = [107, 106, 108]; %phonons to look at, ordered by energy [A2, B2, B1] [5/107, 106, 108]
scale = 2; %Max scale factor in Angstrom sqrt(amu)
gridsize = 9; % Grid size
phase = 0; %phase factor
qpt = 1; %qpt to look at
supercell = [1, 1, 1]; %supercell to construct

%Load relevant quantities

pm = 1822.89; % 1 amu in electron mass units
b2A = 1.88973; % 1 Angstrom in Bohr radii

eigend_abi = ncread(filen, 'phdispl_cart'); %Eigendisplacements - switch to include non-analyticity if necessary
qpoints = ncread(filen, 'qpoints'); %List of qpoints
mass_amu = ncread(filen, 'atomic_mass_units'); %atomic masses in amu
species = ncread(filen, 'atom_species'); %atomic species labels
phfrq = ncread(filen, 'phfreqs'); %Phonon energies in eV
rprim = transpose(ncread(filen, 'primitive_vectors')); % primitive lattice vectors in bohr
xred = transpose(ncread(filen, 'reduced_atom_positions')); %reduced coordinates

natom = size(species,1); %number of atoms
nphonon = size(eigend_abi,2); %number of phonon bands
ncells = supercell(1).*supercell(2).*supercell(3);

mass = mass_amu.*pm; %convert masses to emu

%Reformat eigendisplacements (bohr) and get nornalized ones

for phonon = 1:nphonon
    eigend_abi_complex(:,phonon) = complex(eigend_abi(1,:,phonon,qpt),eigend_abi(2,:,phonon,qpt));
    eigend_abi_normed(:,phonon) = eigend_abi_complex(:,phonon)/norm(eigend_abi_complex(:,phonon));
    for atom = 1:natom
        amf = sqrt((mass(species(atom))));
        for dir = 1:3
            eigendr(atom,dir,phonon) = eigend_abi(1,(atom-1)*3+dir,phonon,qpt); %real part
            eigendi(atom,dir,phonon) = eigend_abi(2,(atom-1)*3+dir,phonon,qpt); %imaginary part
            eigend(atom,dir,phonon) = complex(eigendr(atom,dir,phonon),eigendi(atom,dir,phonon));
            eivecr(atom,dir,phonon) = eigendr(atom,dir,phonon)*amf;
            eiveci(atom,dir,phonon) = eigendi(atom,dir,phonon)*amf;
            eivec(atom,dir,phonon) = complex(eivecr(atom,dir,phonon),eiveci(atom,dir,phonon));
            eigend_normed(atom,dir,phonon) = eigend_abi_normed((atom-1)*3+dir,phonon,qpt);
        end
    end
end

%% Get cartesian coordinates

for atom = 1:natom
    xcart(atom,:) = xred(atom,1)*rprim(1,:)+xred(atom,2)*rprim(2,:)+xred(atom,3)*rprim(3,:);
end

%% Make supercell

rprim_s = transpose(supercell).*rprim;

qvec = supercell.^-1;
qvec(qvec == 1) = 0;

atom_s = 0;

[~, j, k] = deal(1);
Q = linspace(-2+0.000001,scale,gridsize+0.000001);
   
for cellx = 1:supercell(1)
    cellindex(1) = cellx;
    for celly = 1:supercell(2)
        cellindex(2) = celly;
        for cellz = 1:supercell(3)
            cellindex(3) = cellz;
            for atom = 1:natom
                atom_s = atom_s + 1;
                amf = sqrt((mass(species(atom))));
                for dir = 1:3
                    xred_s(atom_s,dir) = xred(atom,dir)/supercell(dir)+(cellindex(dir)-1)/supercell(dir);
                end
                xcart_s(atom_s,:) = xred_s(atom_s,1)*rprim_s(1,:)+xred_s(atom_s,2)*rprim_s(2,:)+xred_s(atom_s,3)*rprim_s(3,:);
                species_s(atom_s) = species(atom);
                for i = 1:gridsize
                    ucart1(atom_s,:) = Q(i)*pm^0.5/b2A/amf*eigend_normed(atom,:,phonon_index(1))*exp(1i*(2*pi*dot(qvec,xred_s(atom_s,:).*transpose(supercell(:)))-phase));
                    for j = 1: gridsize
                        ucart2(atom_s,:) = Q(j)*pm^0.5/b2A/amf*eigend_normed(atom,:,phonon_index(2))*exp(1i*(2*pi*dot(qvec,xred_s(atom_s,:).*transpose(supercell(:)))-phase));
                        for k = 1:gridsize
                            ucart3(atom_s,:) = Q(k)*pm^0.5/b2A/amf*eigend_normed(atom,:,phonon_index(3))*exp(1i*(2*pi*dot(qvec,xred_s(atom_s,:).*transpose(supercell(:)))-phase));
                            xcart_sp(atom_s,:,i,j,k) = xcart_s(atom_s,:) + real(ucart1(atom_s,:)) + real(ucart2(atom_s,:)) + real(ucart3(atom_s,:));
                        end
                    end
                end
            end
        end
    end
end

%% Make some part of the Abinit input

delete mydiary.txt

diary mydiary.txt
diary_folder = pwd; 

format long

dataset = 1;

for i = 1:gridsize
    for j = 1:gridsize
        for k = 1: gridsize
            disp(['xcart', num2str(dataset)])
            disp(xcart_sp(:,:,i,j,k))
            dataset = dataset + 1;
        end
    end
end
dataset = dataset -1;


disp(['ndtset ', num2str(dataset)])
disp(['jdtset '])

nbreak = 25;

for i = 1:floor(dataset/nbreak)
    if i < floor(dataset/nbreak)
        disp(num2str(linspace(1+(i-1)*nbreak,i*nbreak,nbreak)))
    else
        disp(num2str(linspace(1+(i-1)*nbreak,i*nbreak + mod(dataset,nbreak),nbreak + mod(dataset,nbreak))))
    end
end

diary off

