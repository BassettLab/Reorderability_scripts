%% Making all reps for Julia
% Last update 12.17.18


%%%%%% NEED TO ENSURE ALL THREE SECTIONS USING THE SAME GRAPH MODEL %%%%%%

clear

%% Parameters

% filename and params
graphname = 'NF_RG_ep015_1218';


p = 0.4;
d = 3;
ep = 0.15;
a = 1;
maxProb = 0.7;
commSize = 30;
p_in = .6;     %.65 
p_out = 0.4;    %.2
nNodes = 70;
nReps = 1000;
nGraphs = 1;
affinityVec = randperm(nNodes);
alpha = 4;
beta = 1;
maxDeg = 50;
M = 18;     %avg degree is 36 of semantic all
m = 4;
m_0 = 4;



%% Original ordering

filename = graphname;

% Do multiple reps
nReps = 1000;
nGraphs = 1;
adj_array = zeros(nNodes,nNodes,nReps);
jadj_array = adj_array;
s_0_array = zeros(nReps,nNodes);
badj_array = adj_array;
svec_array = cell(nReps,1);
adjn_array = zeros(nNodes,nNodes,nReps);
%load('Stey_unif_f.mat')

for rep = 1:nReps
    
    %aVec = tan(pi*(1:120));
    %affinityVec = randperm(nNodes);
    % Make model
    %[s_0,badj] = NF_ER(nNodes,p);
    %[s_0,badj] = NF_proportionalProb(nNodes,d);
    %[s_0,badj] = NF_fn_polyneg(nNodes,d);
    %[s_0,badj] = NF_fn_frac(nNodes,d);
    [s_0,badj] = NF_fn_abssin(nNodes,[0 2*pi],a);
    %[s_0,badj] = NF_RG(nNodes,ep,d);
    %[s_0,badj] = NF_fn_ln(nNodes);
    %badj = Steyvers_prefAttachUnif(nNodes,M);
    %s_0 = 1:nNodes;
    %badj = badj_array(:,:,rep);
    % s_0 = randperm(nNodes);
    %[s_0,badj] = NF_ER_try2(nNodes,p);
    %[s_0,badj] = NF_randDeg(nNodes,maxDeg);
    %[s_0,badj,commvec] = NF_modular(nNodes,p_in,p_out,commSize);
    %[s_0,badj] = NF_randProb(nNodes,maxProb);
    %[s_0, badj] = NF_affinity(nNodes,affinityVec,a);
    %[s_0,badj,commvec] = NF_modular_affinity(nNodes,p_in,p_out,commSize,affinityVec);
    %load('semanticNorder.mat')
    %[s_0,badj] = NF_modular(nNodes,p_in,p_out,commSize);
    %[s_0] = NF_distanceFromV0_bin(badj,v_0);
    %[s_0] = NF_degree(badj);
    %[ badj, s_0 ] = prefAttachBA( nNodes, m_0,m);
    %[ s_0, badj ] = NF_spatialgrowth(nNodes,beta,alpha);
    
    % Turn nw graph into edge-weighted netowrk
    [jadj,adj,adjn] = makeNodeOrderAdj2(s_0,badj);

    % Storage
    jadj_array(:,:,rep) = jadj;
    adj_array(:,:,rep) = adj;
    s_0_array(rep,:) = s_0;
    badj_array(:,:,rep) = badj;
    %svec_array{rep,1} = svec_out;
    adjn_array(:,:,rep) = adjn;
    
end


save(sprintf('Results/%s_forJul',filename),'jadj_array')
save(sprintf('Results/%s',filename),'badj_array','adj_array','s_0_array',...
    'adjn_array','jadj_array','nReps','nGraphs','nNodes')

disp('Done with original reps :)')





%% Making global reorderings

filename = sprintf('%s_glob',graphname);

% Do multiple reps
nGraphs = 100;
nReps = 100+1;
adj_array = zeros(nNodes,nNodes,nReps*nGraphs);
jadj_array = adj_array;
s_0_array = zeros(nReps*nGraphs,nNodes);
badj_array = adj_array;
svec_array = cell(nReps,1);
adjn_array = zeros(nNodes,nNodes,nReps*nGraphs);
%load('Stey_unif_f.mat')


for graphn = 1:nGraphs
    %aVec = tan(pi*(1:120));
    %affinityVec = randperm(nNodes);
    % Make model
    %[s_0,badj] = NF_proportionalProb(nNodes,d);
    %[s_0,badj] = NF_ER(nNodes,p);
    %[s_0,badj] = NF_fn_polyneg(nNodes,d);
    %[s_0,badj] = NF_fn_frac(nNodes,d);
    %[s_0,badj] = NF_fn_abssin(nNodes,[0 2*pi],a);
    %[s_0,badj] = NF_RG(nNodes,ep,d);
    %badj = Steyvers_prefAttachUnif(nNodes,M);
    %[s_0,badj] = NF_fn_ln(nNodes);
    % s_0 = 1:nNodes;
    %badj = badj_array(:,:,rep);
    % s_0 = randperm(nNodes);
    %[s_0,badj] = NF_ER_try2(nNodes,p);
    %[s_0,badj] = NF_randDeg(nNodes,maxDeg);
    %[s_0,badj,commvec] = NF_modular(nNodes,p_in,p_out,commSize);
    %[s_0,badj] = NF_randProb(nNodes,maxProb);
    %[s_0, badj] = NF_affinity(nNodes,affinityVec,a);
    %[s_0,badj,commvec] = NF_modular_affinity(nNodes,p_in,p_out,commSize,affinityVec);
    %load('semanticNorder.mat')
    %[s_0,badj] = NF_modular(nNodes,p_in,p_out,commSize);
    %[s_0] = NF_distanceFromV0_bin(badj,v_0);
    %[s_0] = NF_degree(badj);
    %[ badj ] = prefAttachUnif( nNodes, M );
    %[ badj ] = modularGrowth_Toivonen( nNodes, N_0, n_sec, p_i );
    %[badj] = Ring_Lattice(nNodes,K);
    %badj(badj>0) = 1;
    %s_0 = 1:nNodes;
    [ badj, s_0 ] = prefAttachBA( nNodes, m_0,m);
    %[ s_0, badj ] = NF_spatialgrowth(nNodes,beta,alpha);
   
    % Turn nw graph into edge-weighted netowrk
    [jadj,adj,adjn] = makeNodeOrderAdj2(s_0,badj);


    % Storage
    jadj_array(:,:,1+(graphn-1)*nReps) = jadj;
    adj_array(:,:,1+(graphn-1)*nReps) = adj;
    s_0_array(1+(graphn-1)*nReps,:) = s_0;
    badj_array(:,:,1+(graphn-1)*nReps) = badj;
    %svec_array{rep,1} = svec_out;
    adjn_array(:,:,1+(graphn-1)*nReps) = adjn;
    
    
 %% Now permute the node order nReps times

  for rep = 2:nReps
 
    
     s_reord = randperm(nNodes);
     
     % Turn nw graph into edge-weighted netowrk
    [jadj,adj,adjn] = makeNodeOrderAdj2(s_reord,badj);
    

    % Storage
    jadj_array(:,:,rep+(graphn-1)*nReps) = jadj;
    adj_array(:,:,rep+(graphn-1)*nReps) = adj;
    s_0_array(rep+(graphn-1)*nReps,:) = s_reord;
    badj_array(:,:,rep+(graphn-1)*nReps) = badj;
    %svec_array{rep,1} = svec_out;
    adjn_array(:,:,rep+(graphn-1)*nReps) = adjn;
 
 
 end

end
save(sprintf('Results/%s_forJul',filename),'jadj_array')
save(sprintf('Results/%s',filename),'badj_array','adj_array','s_0_array',...
    'adjn_array','jadj_array','nGraphs','nReps','nNodes')

fprintf('Saved %s_forJul.mat and %s.mat\n',filename,filename) 
disp('Done with global reorderings :)')


%% Make local reorderings

filename = sprintf('%s_local',graphname);

% Do multiple reps
nGraphs = 20;
nNodes = 70;
nReps = nchoosek(nNodes,2);
adj_array = zeros(nNodes,nNodes,nReps*nGraphs);
jadj_array = adj_array;
s_0_array = zeros(nReps*nGraphs,nNodes);
badj_array = adj_array;
svec_array = cell(nReps,1);
adjn_array = zeros(nNodes,nNodes,nReps*nGraphs);
%load('Stey_unif_f.mat')

iter = 1;
for graphn = 1:nGraphs
    %aVec = tan(pi*(1:120));
    %affinityVec = randperm(nNodes);
    % Make model
    %[s_0,badj] = NF_proportionalProb(nNodes,d);
    %[s_0,badj] = NF_ER(nNodes,p);
    %[s_0,badj] = NF_fn_polyneg(nNodes,d);
    %[s_0,badj] = NF_fn_frac(nNodes,d);
    %[s_0,badj] = NF_fn_abssin(nNodes,[0 2*pi],a);
    [s_0,badj] = NF_RG(nNodes,ep,d);
    %[s_0,badj] = NF_fn_ln(nNodes);
    %badj = Steyvers_prefAttachUnif(nNodes,M);
    %s_0 = 1:nNodes;
    %badj = badj_array(:,:,rep);
    % s_0 = randperm(nNodes);
    %[s_0,badj] = NF_ER_try2(nNodes,p);
    %[s_0,badj] = NF_randDeg(nNodes,maxDeg);
    %[s_0,badj,commvec] = NF_modular(nNodes,p_in,p_out,commSize);
    %[s_0,badj] = NF_randProb(nNodes,maxProb);
    %[s_0, badj] = NF_affinity(nNodes,affinityVec,a);
    %[s_0,badj,commvec] = NF_modular_affinity(nNodes,p_in,p_out,commSize,affinityVec);
    %load('semanticNorder.mat')
    %[s_0,badj] = NF_modular(nNodes,p_in,p_out,commSize);
    %[s_0] = NF_distanceFromV0_bin(badj,v_0);
    %[s_0] = NF_degree(badj);
    %[ badj, s_0 ] = prefAttachBA( nNodes, m_0,m);
    %[ s_0, badj ] = NF_spatialgrowth(nNodes,beta,alpha);
    %[ badj ] = modularGrowth_Toivonen( nNodes, N_0, n_sec, p_i );
    %s_0 = 1:nNodes;
    
   
    % Turn nw graph into edge-weighted netowrk
    [jadj,adj,adjn] = makeNodeOrderAdj2(s_0,badj);


    % Storage
    jadj_array(:,:,iter) = jadj;
    adj_array(:,:,iter) = adj;
    s_0_array(iter,:) = s_0;
    badj_array(:,:,iter) = badj;
    %svec_array{rep,1} = svec_out;
    adjn_array(:,:,iter) = adjn;
    
    iter = iter+1;
    
 % Now permute the node order nReps times
  rep = 2;
  for i1 = 1:nNodes
      for j1 = (i1+1):nNodes
          [ s_reord ] = swapij(s_0,i1,j1);
 
    
          % Turn nw graph into edge-weighted netowrk
          [jadj,adj,adjn] = makeNodeOrderAdj2(s_reord,badj);
          
          % Storage
          jadj_array(:,:,iter) = jadj;
          adj_array(:,:,iter) = adj;
          s_0_array(iter,:) = s_reord;
          badj_array(:,:,iter) = badj;
          %svec_array{rep,1} = svec_out;
          adjn_array(:,:,iter) = adjn;
          
          rep = rep+1;
          iter = iter+1;
      end
  end

end
save(sprintf('Results/%s_forJul',filename),'jadj_array')
save(sprintf('Results/%s',filename),'badj_array','adj_array','s_0_array',...
    'adjn_array','jadj_array','nGraphs','nReps','nNodes')

fprintf('Saved %s_forJul.mat and %s.mat\n',filename,filename) 
disp('Done with local reordering :)')
