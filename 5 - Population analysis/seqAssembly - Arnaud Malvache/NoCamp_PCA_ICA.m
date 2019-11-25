function ICout=NoCamp_PCA_ICA(Img2)

% ********************************************************************** %
% MAIN - PCA ICA FOR CELLS ACTIVITY DETECTION
% ********************************************************************** %

% 1. Load data
% 2. Pre-processing 
% 3. Define parameters for PCA - ICA
% 4. PCA
% 5. ICA
% 6. Plots

% Movie = Img and M
% Img is the 3D matrix heigh X width X time
% M is the 2D matrix pixels X time

% Img1, M1 is the original movie 
% Img2, M2 is the movie after pre-processing (blurr)
% Img3, M3 is the movie after pre-processing (background correction)
% Img4, M4 is the centralized movie for PCA 






%% Load data --------------------------------------------------------------


% -- Extract a new movie -> Img0
% File0='a_C00_T0000.ome.tif';
% Path0='E:\Nathalie\Data\h63_20130702\';
% Img1 = LoadImage2(File0,Path0); 
% Img1 = double(Img1);

% -- Or load a movie -> Img1
% load Film0CorMov
% Img1 = double(Img1);

% -- Or load a blurred movie -> Img2
% load ImgBlurr3


% -- Movie size - number of time frame and number of pixels (height,width)
if exist('Img1','var')
    [Nh0,Nw0,Nt] = size(Img1);
elseif exist('Img2','var')
    [Nh0,Nw0,Nt] = size(Img2);
end

    

    
%% Pre-processing ---------------------------------------------------------


% -- Gaussian blurr
if exist('Img2','var')
    
else
    % Blurr image
    Img2 = GaussBlurXY(Img1,Nh0/3,Nw0/3);
end
% Img2=double(Img2);

% -- Background correction
% matlabpool
% Img3 = CorrBackground(Img2,3,1); % correction using mean image
% Img3 = CorrBackground(Img2,3,0); % correction at each time frame
% matlabpool close
    

% -- Reshape for PCA
if exist('Img3','var')
    M4 = reshape(Img3,Nh0*Nw0,Nt);
else
    M4 = reshape(Img2,Nh0*Nw0,Nt);
end






%% Parameters for PCA - ICA -----------------------------------------------


% -- Auto-cut the principal components
b_cutpc = 1;
% -- Weight eigenvectors with eigenvalues
b_weight = 1;
% -- Number of principal components 
Npc_max = 300;
% -- Number of indepedent components
Nic_max = 300;
% -- Number max of steps in ICA algorithm give you information in space
Nstep_max = 1000;
% -- Coefficient for spatio-temporal ICA (0=temporal, 1=spatial)
mu_ica = 1;






%% Step 1: PCA ------------------------------------------------------------


disp('*** PCA ***')

if Npc_max>min(Nh0*Nw0,Nt)
    disp('Warning: more principal components than the minimal dimension. Change to the minimal dimension.')
    Npc_max = min(Nh0*Nw0,Nt);
end
disp(['PCA with ', num2str(Npc_max),' principal components'])

% tic
[M4,PCout] = step1_PCA(M4,Npc_max);
% toc
disp(['Found ', num2str(length(PCout.pc_EigVal)),' principal components'])   


pc_EigVal = PCout.pc_EigVal;            % eigenvalues
pc_XFilter = PCout.pc_XFilter;          % spatial eigenvectors
pc_TimeCourse = PCout.pc_TimeCourse;    % temporal eigenvectors


if b_cutpc
    
    % -- Auto-cut the PCs -> remove pc with sqrt(eigval) > std noise
    M4std =  mean(std(M4,0,2));
    K_keep = [1 find(sqrt(pc_EigVal)>M4std,1,'last')];
    
    % -- Or choose manually 
    % K_keep = [1 1500]
    
    pc_EigVal = PCout.pc_EigVal(K_keep(1):K_keep(2));
    pc_XFilter = PCout.pc_XFilter(:,K_keep(1):K_keep(2));
    pc_TimeCourse = PCout.pc_TimeCourse(:,K_keep(1):K_keep(2));
    
end

Npc = length(pc_EigVal);
disp(['Keep ', num2str(Npc), ' out of ', num2str(Npc_max), ' eigenvalues'])


if b_weight
     
    disp('Working with weigthed data')
    mat_SVD = diag(sqrt(max(Nh0*Nw0,Nt)*pc_EigVal));
    pc_XFilter = pc_XFilter*(mat_SVD).^(1/2);
    pc_TimeCourse = pc_TimeCourse*(mat_SVD).^(1/2);
    
else % Mwhite = pc_XFilter*pc_TimeCourse';
    disp('Working with non-weigthed data')
end




%% Step 2: ICA ------------------------------------------------------------


disp('*** ICA ***')

if Nic_max>Npc
    disp('Warning: more independent components than principal components. Change to the number of PCs.')
    Nic_max = Npc;
end
disp(['ICA with ', num2str(Nic_max),' independent components'])

% -- Initialize unmixing matrix
UnmixingMat_init = eye(Nic_max,Npc);

% tic 
[ ICout ] = step2_ICA(mu_ica,UnmixingMat_init,pc_XFilter,pc_TimeCourse,Nic_max,Nstep_max,[]);  
% toc