clear
load kiind_fileok/EIds.mat;
load kiind_fileok/NodeIds.mat;
load kiind_fileok/Np.mat;
load kiind_fileok/Tnode.mat;
load kiind_fileok/TP.mat;
load kiind_fileok/VelX.mat;
load kiind_fileok/VelY.mat;
load kiind_fileok/VelZ.mat;

%% Filtering bad nodes
 %Creating Temperature gradient (normalised to v_length)
 
   dt=diff(TP');  % Numerikusan szeretn�m deriv�lni a h�m�rs�kletet az id�vel, ehhez a mennyis�gek k�l�nbs�g�re van sz�ks�gem
   dt(:,1)=[];
    dt(1,:)=[];

   
 
   dT=diff(Tnode);  
   gradient=dT./dt; 
  
   
   [m,n]=size(gradient);
 
 % Setting the target: If any node at any time-step is higher than this
 % value, refining is necessary 
 

  gradient(gradient>10000)=0;
  linear=abs(gradient(:));
  line=sort(linear);
  target=prctile(line, 99);   % Ez a percentilis �rt�k (99) azt mondja meg, hogy melyik a csom�pontok 1%-a, amiben a legnagyobb a differencia.
  
 
  bad_array = abs(gradient)>target;
  bad_row=any(bad_array,1);
  bad_nodes = find(bad_row);
  
  
   
  
%% Defining the bad elements derived from bad nodes
  
 % Handling 'cell' format
 
    maxLengthCell=max(cellfun('size',EIds,2)); 
    for i=1:length(EIds)
        for j=cellfun('size',EIds(i),2)+1:maxLengthCell
             EIds{i}(j)=-1 ;   
        end
    end
   A=cell2mat(EIds);
   Elements=A'+1;
   
   [o,p]=size(Elements);
   
 % Sorting bad elements
 
   bad_Earray=Elements(:,(bad_nodes));
   bad_Erow=bad_Earray(:);
   bad_Elements=unique(bad_Erow);
   bad_Elements(1:1)=[];
   
   
 %% Placing point of balance into bad elements

 % +1 is necessary, for sometimes the IDs start by 0 instead of 1
 
  NodeIds=NodeIds-min(min(NodeIds))+1;
  
  bad_Elements_nodes=NodeIds(:,bad_Elements');
  NpT=Np';
  
 % Centre/Point of balance 
 
  Element_positions=NpT(:,bad_Elements_nodes);

  
  dimension=size(Element_positions,2)/4;
  Element_coordinates=reshape(Element_positions, 3,4,dimension);
  
  PoB=mean(Element_coordinates,2);
  PoB_perfect=reshape(PoB,3,dimension);
%   PoB_perfect(3,:)=PoB_perfect(3,:)+0.00002+(0.00007-0.00002)*rand(1,dimension);
%   PoB_perfect(2,:)=PoB_perfect(2,:)+0.00002+(0.00007-0.00002)*rand(1,dimension);
%   PoB_perfect(1,:)=PoB_perfect(1,:)+0.00002+(0.00007-0.00002)*rand(1,dimension);
  
   %% The new "cloud of nodes"
   
   % Handling 'cell' format
 
    maxLengthCell=max(cellfun('size',EIds,2)); 
    for i=1:length(EIds)
        for j=cellfun('size',EIds(i),2)+1:maxLengthCell
             EIds{i}(j)=-1 ;   
        end
    end
   A=cell2mat(EIds);
   Elements=A'+1;
   
   [o,p]=size(Elements);
  
  Points=[NpT PoB_perfect];
  NodeIds(:, bad_Elements)=[];
  nn=size(PoB_perfect);
  baddies=nn(2);
  newNodes=(p+1):1:(p+baddies);
  bad_cd=1:1:baddies;
  
  newElements=[bad_Elements_nodes(2, bad_cd) bad_Elements_nodes(1, bad_cd)  bad_Elements_nodes(1, bad_cd) bad_Elements_nodes(2, bad_cd);
                bad_Elements_nodes(1, bad_cd) bad_Elements_nodes(2, bad_cd) bad_Elements_nodes(4, bad_cd) bad_Elements_nodes(3, bad_cd);
                bad_Elements_nodes(3, bad_cd) bad_Elements_nodes(4, bad_cd) bad_Elements_nodes(3, bad_cd) bad_Elements_nodes(4, bad_cd);
                newNodes newNodes newNodes newNodes];
 % new_Elements=... "var�zsl�s": Az �j elemek defini�l�sa nem annyira trivi�lis. A v�geselemes szoftverben (Moldflownak) az elemeket adott �gynevezett k�r�lj�r�si ir�nynak kell megfeleltetni. Ez azt jelenti, hogy p�ld�ul az 1,2,3,4-es sz�m� csom�pontok, amikor defini�ljuk vel�k az elemeket nem mindegy, hogy milyen sorrendben vannak felt�ntetve.           
  ConnectivityList=[NodeIds newElements];
  
%% Figures

   position=Np(bad_nodes,:);
  
   X=Np(:,1);
   Y=Np(:,2);
   Z=Np(:,3);
     
   Xp=position(:,1);
   Yp=position(:,2);
   Zp=position(:,3);
   
   Xc=PoB_perfect(1,:);
   Yc=PoB_perfect(2,:);
   Zc=PoB_perfect(3,:);
   
   hold on
   scatter3(X,Y,Z, 'b')
   scatter3(Xp,Yp,Zp, 'r')
   scatter3(Xc,Yc,Zc, 'g')
  

%% �rk�l�s UDM-be % Erre nem lesz egyel�re sz�ks�getek.
  
%  a=1:size(ConnectivityList, 2);
%    
% formatSpec0='\n NOT4{%d}}\n \n';
% p=size(Points, 2);
% g=1:1:p;
% str0=sprintf(formatSpec0, p);
% 
% formatSpec = 'NODE{%u 0 2 1 0.1  %0.6e  %0.6e  %0.6e}\n';
% formatSpec2 = 'TET4{%u 0 4 27  0 "" 50400 1       %u       %u       %u       %u} \n';
%    
% fid=fopen('egyelem.udm');
% % save('finom.udm');
% % fclose(fi);
% fid=fopen('egyelem.udm', 'a+t');
% fprintf(fid, str0);
% fprintf(fid, formatSpec,[g; Points]);
% fprintf(fid, formatSpec2,[a; ConnectivityList]);
% fprintf(fid, 'ENDF');
% fclose(fid);   

  
  
 
    
      
   

      
  

  
 
           
 
 
  
 
  
  


  

      
 
 
  