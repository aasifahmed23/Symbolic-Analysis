
function [V,X,SymString] = scam_v5(fid,Level)

C=textscan(fid,'%s %d16 %d16 %f32 %d16 %s %s %s','EmptyValue', -1); 
Name=C{1};
N1=C{2};
N2=C{3};
arg3=C{4};
N4=C{5};
Type=C{6};
L=C{7};
W=C{8};
[rows cols]=size(C{1});
%
for i=1:rows
if ( ~isempty(L{i}))
    L_new=textscan(L{i},'L=%f');
    L(i)=L_new;
end

if ( ~isempty(W{i}))
    W_new=textscan(W{i},'W=%f');
    W(i)=W_new;
end
end


tic
% Initialize
numElem=0;  %Number of passive elements.
numV=0;     %Number of independent voltage sources
numO=0;     %Number of op amps
numI=0;     %Number of independent current sources
numNode=0;  %Number of nodes, not including ground (node 0).
numM=0;     %Number of MOS devices

% Parse the input file
for i=1:length(Name),
    switch(Name{i}(1)),
        case {'R','L','C'},
            numElem=numElem+1;
            Element(numElem).Name=Name(i);
            Element(numElem).Node1=N1(i);
            Element(numElem).Node2=N2(i);
            try
                Element(numElem).Value=arg3(i);
            catch
                Element(numElem).Value=nan;
            end
        case 'V',
            numV=numV+1;
            Vsource(numV).Name=Name(i);
            Vsource(numV).Node1=N1(i);
            Vsource(numV).Node2=N2(i);
            try
                Vsource(numV).Value=arg3(i);
            catch
                Vsource(numV).Value=nan;
            end
        case 'O',
            numO=numO+1;
            Opamp(numO).Name=Name(i);
            Opamp(numO).Node1=N1(i);
            Opamp(numO).Node2=N2(i);
            Opamp(numO).Node3=arg3(i);
        case 'I'
            numI=numI+1;
            Isource(numI).Name=Name(i);
            Isource(numI).Node1=N1(i);
            Isource(numI).Node2=N2(i);
            try
                Isource(numI).Value=arg3(i);
            catch
                Isource(numI).Value=nan;
            end
        case 'M'
            numM=numM+1;
            MOS(numM).Name=Name(i);
            MOS(numM).Node1=N1(i);
            MOS(numM).Node2=N2(i);
            MOS(numM).Node3=arg3(i);
            MOS(numM).Node4=N4(i);
            MOS(numM).Type=Type(i);
            MOS(numM).W=W(i);
            MOS(numM).L=L(i);
    end
    numNode=max(N1(i),max(N2(i),numNode));
end

% Preallocate all of the cell arrays #################################
G=cell(numNode,numNode);
V=cell(numNode,1);
I=cell(numNode,1);
if ((numV+numO)~=0),
    B=cell(numNode,numV+numO);
    C=cell(numV+numO,numNode);
    D=cell(numV+numO,numV+numO);
    E=cell(numV+numO,1);
    J=cell(numV+numO,1);
end
%Done preallocating cell arrays -------------------------------------

% Fill the G matrix ##################################################
%Initially, make the G Matrix all zeros.
[G{:}]=deal('0');

%Now fill the G matrix with conductances from netlist
for i=1:numElem,
    n1=Element(i).Node1;
    n2=Element(i).Node2;
    %Make up a string with the conductance of current element.
    switch(Element(i).Name{1}(1)),
        case 'R',
            g = ['1/' Element(i).Name{1}];
        case 'L',
            g = ['1/s/' Element(i).Name{1}];
        case 'C',
            g = ['s*' Element(i).Name{1}];
    end
    
    %If neither side of the element is connected to ground
    %then subtract it from appropriate location in matrix.
    if (n1~=0) && (n2~=0),
        G{n1,n2}=[ G{n1,n2} '-' g];
        G{n2,n1}=[ G{n2,n1} '-' g];
    end
    
    %If node 1 is connected to graound, add element to diagonal
    %of matrix.
    if (n1~=0),
        G{n1,n1}=[ G{n1,n1} '+' g];
    end
    %Ditto for node 2.
    if (n2~=0),
        G{n2,n2}=[ G{n2,n2} '+' g];
    end
    
    %Go to next element.
    %     i=i+4;
end
F=cell(4,4);
[F{:}]=deal('0');
F1=cell(4,4);
[F1{:}]=deal('0');
% MNA Stamping of MOSFETS

for i=1:numM, % for each MOSFETS
    
    n1=MOS(i).Node1; % find the drain node
    n2=MOS(i).Node2; % find the gate node
    n3=MOS(i).Node3; % find the source node
    n4=MOS(i).Node4; % find the substrate node
    n=[n1 n2 n3 n4]; 
    % n=[3 3 2 1]; 
    % Matlab Warning: Concatenation with dominant (left-most) integer class 
    % may overflow other operands on conversion to return class.
    % (Type "warning off MATLAB:concatenation:integerInteraction" to suppress this warning.)

%1#####calculate the device stamp 'F' for un-Shorted,un-Grounded nodes#####
% Level=3; % Switch for choosing the level of MOSFET model
switch Level
    case 4 % MOSFET model contains gm,gmb,gds,Cgd,Cgs,Cdb,Csb,Cgb (8 small signal parameters)
   g11=['s*' 'Cgd' '_' MOS(i).Name{1} '+' 's*' 'Cdb' '_' MOS(i).Name{1} '+' 'gds' '_' MOS(i).Name{1}];
   g22=['s*' 'Cgd' '_' MOS(i).Name{1} '+' 's*' 'Cgs' '_' MOS(i).Name{1} '+' 's*' 'Cgb' '_' MOS(i).Name{1}];
   g33=['s*' 'Cgs' '_' MOS(i).Name{1} '+' 's*' 'Csb' '_' MOS(i).Name{1} '+' 'gm' '_' MOS(i).Name{1} '+' 'gds' '_' MOS(i).Name{1} '+' 'gmb' '_' MOS(i).Name{1}];
   g44=['s*' 'Cdb' '_' MOS(i).Name{1} '+' 's*' 'Csb' '_' MOS(i).Name{1} '+' 's*' 'Cgb' '_' MOS(i).Name{1}];
   F{1,1}=g11;
   F{1,2}=['gm' '_' MOS(i).Name{1} '-' 's*' 'Cgd' '_' MOS(i).Name{1}];
   F{1,3}=['-' 'gds' '_' MOS(i).Name{1} '-' 'gm' '_' MOS(i).Name{1} '-' 'gmb' '_' MOS(i).Name{1}];
   F{1,4}=['gmb' '_' MOS(i).Name{1} '-' 's*' 'Cdb' '_' MOS(i).Name{1}];
   F{2,1}=['-' 's*' 'Cgd' '_' MOS(i).Name{1}];
   F{2,2}=g22;
   F{2,3}=['-' 's*' 'Cgs' '_' MOS(i).Name{1}];
   F{2,4}=['-' 's*' 'Cgb' '_' MOS(i).Name{1}];
   F{3,1}=['-' 'gds' '_' MOS(i).Name{1}];
   F{3,2}=['-' 's*' 'Cgs' '_' MOS(i).Name{1} '-' 'gm' '_' MOS(i).Name{1}];
   F{3,3}=g33;
   F{3,4}=['-' 's*' 'Csb' '_' MOS(i).Name{1} '-' 'gmb' '_' MOS(i).Name{1}];
   F{4,1}=['-' 's*' 'Cdb' '_' MOS(i).Name{1}];
   F{4,2}=['-' 's*' 'Cgb' '_' MOS(i).Name{1}];
   F{4,3}=['-' 's*' 'Csb' '_' MOS(i).Name{1}];
   F{4,4}=g44;          
    case 3 % MOSFET model contains gm,gmb,gds,Cgd,Cgs,Cdb,Csb (7 small signal parameters--Cgb is absent)
   g11=['s*' 'Cgd' '_' MOS(i).Name{1} '+' 's*' 'Cdb' '_' MOS(i).Name{1} '+' 'gds' '_' MOS(i).Name{1}];
   g22=['s*' 'Cgd' '_' MOS(i).Name{1} '+' 's*' 'Cgs' '_' MOS(i).Name{1}];
   g33=['s*' 'Cgs' '_' MOS(i).Name{1} '+' 's*' 'Csb' '_' MOS(i).Name{1} '+' 'gm' '_' MOS(i).Name{1} '+' 'gds' '_' MOS(i).Name{1} '+' 'gmb' '_' MOS(i).Name{1}];
   g44=['s*' 'Cdb' '_' MOS(i).Name{1} '+' 's*' 'Csb' '_' MOS(i).Name{1}];
   F{1,1}=g11;
   F{1,2}=['gm' '_' MOS(i).Name{1} '-' 's*' 'Cgd' '_' MOS(i).Name{1}];
   F{1,3}=['-' 'gds' '_' MOS(i).Name{1} '-' 'gm' '_' MOS(i).Name{1} '-' 'gmb' '_' MOS(i).Name{1}];
   F{1,4}=['gmb' '_' MOS(i).Name{1} '-' 's*' 'Cdb' '_' MOS(i).Name{1}];
   F{2,1}=['-' 's*' 'Cgd' '_' MOS(i).Name{1}];
   F{2,2}=g22;
   F{2,3}=['-' 's*' 'Cgs' '_' MOS(i).Name{1}];
   F{2,4}='0';
   F{3,1}=['-' 'gds' '_' MOS(i).Name{1}];
   F{3,2}=['-' 's*' 'Cgs' '_' MOS(i).Name{1} '-' 'gm' '_' MOS(i).Name{1}];
   F{3,3}=g33;
   F{3,4}=['-' 's*' 'Csb' '_' MOS(i).Name{1} '-' 'gmb' '_' MOS(i).Name{1}];
   F{4,1}=['-' 's*' 'Cdb' '_' MOS(i).Name{1}];
   F{4,2}='0';
   F{4,3}=['-' 's*' 'Csb' '_' MOS(i).Name{1}];
   F{4,4}=g44;
    case 2 % MOSFET model contains gm,gmb,gds(3 small signal parameters--Cgd,Cgs,Cdb,Csb,,Cgb are absent)
   g11=['gds' '_' MOS(i).Name{1}];
   g22='0';
   g33=['gm' '_' MOS(i).Name{1} '+' 'gds' '_' MOS(i).Name{1} '+' 'gmb' '_' MOS(i).Name{1}];
   g44='0';
   F{1,1}=g11;
   F{1,2}=['gm' '_' MOS(i).Name{1}];
   F{1,3}=['-' 'gds' '_' MOS(i).Name{1} '-' 'gm' '_' MOS(i).Name{1} '-' 'gmb' '_' MOS(i).Name{1}];
   F{1,4}=['gmb' '_' MOS(i).Name{1}];
   F{2,1}='0';
   F{2,2}=g22;
   F{2,3}='0';
   F{2,4}='0';
   F{3,1}=['-' 'gds' '_' MOS(i).Name{1}];
   F{3,2}=['-' 'gm' '_' MOS(i).Name{1}];
   F{3,3}=g33;
   F{3,4}=['-' 'gmb' '_' MOS(i).Name{1}];
   F{4,1}='0';
   F{4,2}='0';
   F{4,3}='0';
   F{4,4}=g44;
    case 1 % MOSFET model contains gm,gds(3 small signal parameters--gmb,,Cgd,Cgs,Cdb,Csb,,Cgb are absent)
   g11=['gds' '_' MOS(i).Name{1}];
   g22='0';
   g33=['gm' '_' MOS(i).Name{1} '+' 'gds' '_' MOS(i).Name{1}];
   g44='0';
   F{1,1}=g11;
   F{1,2}=['gm' '_' MOS(i).Name{1}];
   F{1,3}=['-' 'gds' '_' MOS(i).Name{1} '-' 'gm' '_' MOS(i).Name{1}];
   F{1,4}='0';
   F{2,1}='0';
   F{2,2}=g22;
   F{2,3}='0';
   F{2,4}='0';
   F{3,1}=['-' 'gds' '_' MOS(i).Name{1}];
   F{3,2}=['-' 'gm' '_' MOS(i).Name{1}];
   F{3,3}=g33;
   F{3,4}='0';
   F{4,1}='0';
   F{4,2}='0';
   F{4,3}='0';
   F{4,4}=g44;
    otherwise % MOSFET model contains gm,gds(3 small signal parameters--gmb,,Cgd,Cgs,Cdb,Csb,,Cgb are absent)  
   g11=['gds' '_' MOS(i).Name{1}];
   g22='0';
   g33=['gm' '_' MOS(i).Name{1} '+' 'gds' '_' MOS(i).Name{1}];
   g44='0';
   F{1,1}=g11;
   F{1,2}=['gm' '_' MOS(i).Name{1}];
   F{1,3}=['-' 'gds' '_' MOS(i).Name{1} '-' 'gm' '_' MOS(i).Name{1}];
   F{1,4}='0';
   F{2,1}='0';
   F{2,2}=g22;
   F{2,3}='0';
   F{2,4}='0';
   F{3,1}=['-' 'gds' '_' MOS(i).Name{1}];
   F{3,2}=['-' 'gm' '_' MOS(i).Name{1}];
   F{3,3}=g33;
   F{3,4}='0';
   F{4,1}='0';
   F{4,2}='0';
   F{4,3}='0';
   F{4,4}=g44; 
        
end % switch ends here
   
     
%2######Finding number of grounded nodes and its position in nodal array 'n'
z=zeros(1,4);
r_grounded=4; % row size initialization
n1=n; %temporary store original node list 'n' in 'n1' for internal processing 
F1=F; %temporary store 'F' in 'F1' for internal processing 
gnd_flag=0;% reset gnd_flag
for j=1:4,
    if (n(j)==0),
        z(j)=1; % captures the grounded nodes
        gnd_flag=1;
        r_grounded=r_grounded-1; % decrement row size for each grounded node
        for m=0:3-j
            n1(j+m)=n(j+1+m); %%possible bug --
        end
     end
end 
%3######Modifying the stamp 'F1' for grounded nodes by removing the 
%3######contribution of grounded nodes from original stamp 'F1'
for j=1:4 % for each node j
if (z(j)==1) % for grounded node case    
  %%%%%%%%%%%%%% LOGIC FOR GROUNDED NODES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for l=0:3-j
        for m=1:4
            F1{(j+l),m}=F1{(j+l+1),m}; % phase 1--row compaction
            F1{(j+l+1),m}=' ';
        end
    end   
    
    for l=0:3-j
        for m=1:4
            F1{m,(j+l)}=F1{m,(j+l+1)}; %phase2--column compaction
            F1{m,(j+l+1)}=' ';
        end
    end
    
    %phase3--storing device stamp for grounded node in 'F2'
    F2=cell(r_grounded,r_grounded);
    for p=1:r_grounded
        for q=1:r_grounded
            F2{p,q}=F1{p,q};
        end
    end    
%%%%%%%%%%%%%%%%%% LOGIC FOR GROUNDED NODES ENDS %%%%%%%%%%%%%%%%%%%%%%%%%%   
end
end

%%4###### node list for grounded nodes
n_grounded=zeros(1,r_grounded);
for j=1:r_grounded
    n_grounded(j)=n1(j);
end

%%5##############Modifying the stamp F1 for shorted nodes
r_shorted=4; %row size initializations
shorts=zeros(3,2);
short_flag=0;
if (gnd_flag==1)
    n=n_grounded;
    F=F2;
    r_shorted=r_shorted-1;
end
n1=n; %temporary store original node list 'n' in 'n1' for internal processing 
F1=F; %temporary store 'F' in 'F1' for internal processing  

for j=1:r_grounded        
        for k=(j+1):r_grounded
            if (n(j)==n(k))&&(n(j)~=0)&&(n(k)~=0) % for shorted non-grounded node case
                r_shorted=r_shorted-1; % decrement row size for each shorted node
                shorts(j,:)=[j k];     % captures the shorted nodes
                short_flag=1;
                for m=0:r_grounded-1-j
                n1(j+m)=n(j+1+m);
                end
%%%%%%%%%%%%%%%%%% LOGIC FOR SHORTED NODES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % phase 1--row merging
           for l=1:r_grounded
               if  (n(l)~=0)
                  if (F{max(k,j),l}(1)~='-')
                  F1{min(k,j),l}=[F{min(k,j),l} '+' F{max(k,j),l}]; % row merging
                  F1{max(k,j),l}='0';
                  else
                  F1{min(k,j),l}=[F{min(k,j),l} F{max(k,j),l}]; % row merging
                  F1{max(k,j),l}='0';
                  end
               end
           end
           %phase 2---column merging 
           for l=1:r_grounded
               if (n(l)~=0)
                   if (F{l,max(k,j)}(1)~='-')
                   F1{l,min(k,j)}=[F{l,min(k,j)} '+' F{l,max(k,j)}]; % column merging 
                   F1{l,max(k,j)}='0';
                   else
                   F1{l,min(k,j)}=[F{l,min(k,j)} F{l,max(k,j)}]; % column merging
                   F1{l,max(k,j)}='0';
                  end
               end
           end
           %phase 3--corner merging
           if (F1{max(k,j),min(k,j)}(1)~='-')
           F1{min(k,j),min(k,j)}=[F1{min(k,j),min(k,j)} '+' F1{max(k,j),min(k,j)}];
           else
               F1{min(k,j),min(k,j)}=[F1{min(k,j),min(k,j)} F1{max(k,j),min(k,j)}];
           end
           F1{max(k,j),min(k,j)}='0';
           
           % phase 4--row compaction phase
           for l=0:r_grounded-1-max(k,j)
             for m=1:r_grounded
               F1{max(k,j)+l,m}=F1{max(k,j)+l+1,m};  %row compaction
               F1{max(k,j)+l+1,m}=' ';
             end
           end
          % phase 5--column compaction phase
          for l=0:r_grounded-1-max(k,j)
           for m=1:r_grounded
            F1{m,max(k,j)+l}=F1{m,max(k,j)+l+1}; %column compaction
            F1{m,max(k,j)+l+1}=' ';
           end
          end 
         %phase 6--storing device stamp for shorted node in 'F2'
         F2=cell(r_shorted,r_shorted);
         for p=1:r_shorted
           for q=1:r_shorted
            F2{p,q}=F1{p,q};
           end
         end
         
        % phase 7--node list for shorted nodes
        n_shorted=zeros(1,r_shorted);
        for p=1:r_shorted
          n_shorted(p)=n1(p);
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%% LOGIC FOR SHORTED NODES ENDS %%%%%%%%%%%%%%%%%%
            end
        end
end

%play ground for final stamping for MOS(i)
%we have actual stamps in 'F'.
%final stamps are in 'F2'.('F1' STAMP IS A SCRATCH PAD ONLY--IT'S USELESS)
%nodelist for grounded nodes are in n_grounded
%nodelist for shorted nodes are in n_shorted
%final nodelist and stamps are in 'n_shorted' & 'F2' of dimension r_shorted
if (gnd_flag==1)
    if(short_flag==1)
        n_final=n_shorted;
        r_final=r_shorted;
    else
        n_final=n_grounded;
        r_final=r_grounded;
    end
else
    if(short_flag==1)
        n_final=n_shorted;
        r_final=r_shorted;
    else
        n_final=n;
        r_final=4;        
    end
end
for k=1:r_final
  for j=1:r_final
      if (F2{k,j}(1)~='-')
    G{n_final(k),n_final(j)}=[G{n_final(k),n_final(j)} '+' F2{k,j}];
      else
    G{n_final(k),n_final(j)}=[G{n_final(k),n_final(j)} F2{k,j}]; 
      end

  end
end

 end
%The G matrix is finished -------------------------------------------

% Fill the I matrix ##################################################
[I{:}]=deal('0');
for j=1:numNode,
    for i=1:numI,
        if (Isource(i).Node1==j),
            I{j}=[I{j} '-' Isource(i).Name{1}];
        elseif (Isource(i).Node2==j),
            I{j}=[I{j} '+' Isource(i).Name{1}];
        end
    end
end
%The I matrix is done -----------------------------------------------

% Fill the V matrix ##################################################
for i=1:numNode,
    V{i}=['v_' num2str(i)];
end
%The V matrix is finished -------------------------------------------
% Fill the B matrix ##################################################
if ((numV+numO)~=0)
        %Initially, fill with zeros.
    [B{:}]=deal('0');
    
    %First handle the case of the independent voltage sources.
    for i=1:numV,           %Go through each independent source.
        for j=1:numNode     %Go through each node.
            if (Vsource(i).Node1==j),       %If node is first node,
                B{j,i}='1';                 %then put '1' in the matrices.
            elseif (Vsource(i).Node2==j),   %If second node, put -1.
                B{j,i}='-1';
            end
        end
    end
    
    %Now handle the case of the Op Amp
    for i=1:numO,
        for j=1:numNode
            if (Opamp(i).Node3==j),
                B{j,i+numV}='1';
            else
                B{j,i+numV}='0';
            end
        end
    end
    %The B matrix is finished -------------------------------------------
   
    
    %%Fill the C matrix ##################################################
    %Initially, fill with zeros.
    [C{:}]=deal('0');
    
    %First handle the case of the independent voltage sources.
    for i=1:numV,           %Go through each independent source.
        for j=1:numNode     %Go through each node.
            if (Vsource(i).Node1==j),       %If node is first node,
                C{i,j}='1';                 %then put '1' in the matrices.
            elseif (Vsource(i).Node2==j),   %If second node, put -1.
                C{i,j}='-1';
            end
        end
    end
    
    %Now handle the case of the Op Amp
    for i=1:numO,
        for j=1:numNode
            if (Opamp(i).Node1==j),
                C{i+numV,j}='1';
            elseif (Opamp(i).Node2==j),
                C{i+numV,j}='-1';
            else
                C{i+numV,j}='0';
            end
        end
    end
    %The C matrix is finished ------------------------------------------
    
    
    %%Fill the D matrix ##################################################
    %The D matrix is non-zero only for CCVS and VCVS (not included
    %in this simple implementation of SPICE)
    [D{:}]=deal('0');
    %The D matrix is finished -------------------------------------------
    
    %Fill the E matrix ##################################################
    %Start with all zeros
    [E{:}]=deal('0');
    for i=1:numV,
        E{i}=Vsource(i).Name{1};
    end
    %The E matrix is finished -------------------------------------------
    
    % Fill the J matrix ##################################################
    for i=1:numV,
        J{i}=['I_' Vsource(i).Name{1}];
    end
    for i=1:numO,
        J{i+numV}=['I_' Opamp(i).Name{1}];
    end
    %The J matrix is finished -------------------------------------------
end  %if ((numV+numO)~=0)

% Form the A, X, and Z matrices (As cell arrays of strings).
if ((numV+numO)~=0),
    Acell=[deal(G) deal(B); deal(C) deal(D)];
    Xcell=[deal(V); deal(J)];
    Zcell=[deal(I); deal(E)];
else
    Acell=[deal(G)];
    Xcell=[deal(V)];
    Zcell=[deal(I)];
end


% Declare symbolic variables #########################################
%This next section declares all variables used as symbolic variables.
%Make "s" a symbolic variable
SymString='syms s ';

%Add each of the passive elements to the list of symbolic variables.
for i=1:numElem,
    SymString=[SymString Element(i).Name{1} ' '];
end

%Add each element of matrix J and E to the list of symbolic variables.
for i=1:numV,
    SymString=[SymString J{i} ' '];
    SymString=[SymString E{i} ' '];
end

%Add each opamp output to the list of symbolic variables.
for i=1:numO,
    SymString=[SymString J{i+numV} ' '];
end

%Add independent current sources to the list of symbolic variables.
for i=1:numI,
    SymString=[SymString Isource(i).Name{1} ' '];
end

%Add independent voltage sources to list of symbolic variables.
for i=1:numNode,
    SymString=[SymString V{i} ' '];
end
for i=1:numM,
    Cgd{i}=['Cgd_' MOS(i).Name{1}];
    Cgs{i}=['Cgs_' MOS(i).Name{1}];
    Cgb{i}=['Cgb_' MOS(i).Name{1}];
    Cdb{i}=['Cdb_' MOS(i).Name{1}];
    Csb{i}=['Csb_' MOS(i).Name{1}];
    gm{i}=['gm_' MOS(i).Name{1}];
    gmb{i}=['gmb_' MOS(i).Name{1}];
    gds{i}=['gds_' MOS(i).Name{1}];
end
%Add MOSFET small signal parameters to list of symbolic variables.
for i=1:numM,
    SymString=[SymString Cgd{i} ' ' Cgs{i} ' ' Cgb{i} ' ' Cdb{i} ' ' Csb{i} ' ' gm{i} ' ' gmb{i} ' ' gds{i} ' '];
end

% Evaluate the string with symbolic variables
eval(SymString);
%Done declaring symbolic variables ----------------------------------

%Create the variables A, X, and Z ###################################
%Right now the matrices Acell, Xcell and Zcell hold cell arrays of 
%strings.  These must be converted to a symbolic array.  This is
%accompplished by creating strings that represent the assignment of
%the symbolic arrays, and then evaluating these strings.

% Create assignments for three arrays
Astring='A=[';
Xstring='X=[';
Zstring='Z=[';

for i=1:length(Acell),     %for each row in the arrays.
    for j=1:length(Acell),      %for each column in matrix A.
        Astring=[Astring ' ' Acell{i,j}]; %Get element from Acell
    end
    Astring=[Astring ';'];          %Mark end of row with semicolon
    Xstring=[Xstring  Xcell{i} ';'];    %Enter element into array X;
    Zstring=[Zstring  Zcell{i} ';'];    %Enter element into array Z;
end
Astring=[Astring '];'];  %Close array assignment.
Xstring=[Xstring '];'];
Zstring=[Zstring '];'];

% Evaluate strings with array assignments.
eval([Astring ' ' Xstring ' ' Zstring])
%Done creating the variables A, X, and Z ----------------------------


% Solve matrrix equation - this is the meat of the algorithm.
V=simplify(inv(A)*Z);

%Evaluate each of the unknowns in the matrix X.
for i=1:length(V),
    eval([char(X(i)) '=' char(V(i)) ';']);
end

disp(sprintf('Done! Elapsed time = %g seconds.\n',toc));
beep;
disp('Solved variables:');
disp(X)


