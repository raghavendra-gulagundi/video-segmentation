function frNO=vsegment(vid_samp,method,th_h,th_l)

g=aviinfo(vid_samp);        % take video information
n=g.NumFrames ;              % get the total number of frames in the video
f=mmreader(vid_samp);       % Create object file for acessing frames

i=1;                        
step=5;      


%*****    START of DCT Mean    *****

if(method==1)
for j=1:step:n              % loop for reading after every "step" number of frames
rf=read(f,j);               % take the jth frame             
r=rf(:,:,1);                % extract the red colour content from frame in 'rf'
g = rf (:,:,2);             % extract the green colour content
b = rf (:,:,3);             % extract the blue colour content
   
[R e]=reshaped(r);          % converts 'r' matrix elements into single row elements and 
G=reshaped(g);              %  'e' is number of zero paded elements during reshaping
B=reshaped(b);
    
MR=mean1(R,e,1);         
MG=mean1(G,e,1);         
MB=mean1(B,e,1);         

MR64(i,:)=MR;               % store the MR,MG,MB of singlw row in the matrix 
MG64(i,:)=MG;
MB64(i,:)=MB;
i=i+1 ;                     % increment the row count of MR64,MG64,MB64 

end
 
Mt(:,1)=(MR64(:,2)+MG64(:,2)+MB64(:,2)+MR64(:,9)+MG64(:,9)+MB64(:,9))/6;   % for mean take both 2nd and 9th AC coefficient
se=size(Mt,1);
Mt=abs(Mt);
gh=2;
for d=1:1:(se-1)
    MS(d,:)=abs(Mt(gh,:)-Mt(d,:));     % taking difference i.e, Current value-Previous value 
    gh=gh+1;
end

%*****    END of DCT Mean    *****



%*****    START of DCT Standard deviation    *****

elseif(method==2) 
for j=1:step:n              
rf=read(f,j);               % take the jth frame 
r=rf(:,:,1);                % extract the red colour content from frame in 'rf'
g = rf (:,:,2);             % extract the green colour content
b = rf (:,:,3);             % extract the blue colour content
   
[R e]=reshaped(r);          % converts 'r' matrix elements into single row elements and 
G=reshaped(g);              %  'e' is numer of zero paded elements lly with g,b
B=reshaped(b);
    
MR=sqrt(nce1(R,e));         % gives the standard deviation of R,G,B column-wise
MG=sqrt(nce1(G,e));         % nce1 gives variance  
MB=sqrt(nce1(B,e));         

MR64(i,:)=MR;               % store the MR,MG,MB of singlw row in the matrix 
MG64(i,:)=MG;
MB64(i,:)=MB;
i=i+1  ;                    % increment the row count of MR64,MG64,MB64 
end
  
MC(:,1) = MR64(:,2) + MG64(:,2) + MB64(:,2);    % take the addition of corresponding 2nd column values

   sz=size(MC,1);           
  op=2;                     
   sz=sz-1;                 
   for h=1:1:sz             
     MS(h,1)= abs(MC(op,1)-MC(h,1)); % take difference of present and previous value in column matrix MC
     op=op+1;                       
   end


%*****    END of DCT Standard deviation     *****



%*****    START of DFT Mean    *****
   
 
  elseif(method==3)
      for j=1:step:n              
rf=read(f,j);               % take the jth frame             
r=rf(:,:,1);                % extract the red colour content from frame in 'rf'
g = rf (:,:,2);             % extract the green colour content
b = rf (:,:,3);             % extract the blue colour content
   
[R e]=reshaped_dft_mg(r);          % converts 'r' matrix elements into single row elements and 
G=reshaped_dft_mg(g);              %  'e' is numer of zero paded elements 
B=reshaped_dft_mg(b);
    
MR=mean1(R,e,1);         % gives the standard deviation of R,G,B column-wise
MG=mean1(G,e,1);         % nce1 is variance function 
MB=mean1(B,e,1);         

MR64(i,:)=MR;               % store the MR,MG,MB of singlw row in the matrix 
MG64(i,:)=MG;
MB64(i,:)=MB;
i=i+1 ;                     % increment the row count of MR64,MG64,MB64 

end
 
Mt(:,1)=(MR64(:,2)+MG64(:,2)+MB64(:,2)+MR64(:,9)+MG64(:,9)+MB64(:,9))/6;
se=size(Mt,1);
Mt=abs(Mt);
gh=2;
for d=1:1:(se-1)
    MS(d,:)=abs(Mt(gh,:)-Mt(d,:));
    gh=gh+1;
end



%*****    END of DFT Mean    *****



%*****    START of DFT Standard deviation    *****
   

      
else
for j=1:step:n              % loop for reading after every "step" number of frames
rf=read(f,j);               % take the jth frame 
               
r=rf(:,:,1);                % extract the red colour content from frame in 'rf'
g = rf (:,:,2);             % extract the green colour content
b = rf (:,:,3);             % extract the blue colour content
   
[R e]=reshaped_dft_mg(r);          % converts 'r' matrix elements into single row elements and 
G=reshaped_dft_mg(g);              %  'e' is numer of zero paded elements 
B=reshaped_dft_mg(b);
    
MR=sqrt(nce1(R,e));         % gives the standard deviation of R,G,B column-wise
MG=sqrt(nce1(G,e));         % nce1 is variance function 
MB=sqrt(nce1(B,e));         

MR64(i,:)=MR;               % store the MR,MG,MB of singlw row in the matrix 
MG64(i,:)=MG;
MB64(i,:)=MB;
i=i+1  ;                    % increment the row count of MR64,MG64,MB64 

end
  
MC(:,1) = MR64(:,2) + MG64(:,2) + MB64(:,2);    % take the addition of corresponding 2nd column values

   sz=size(MC,1);          
  op=2;                    
   sz=sz-1;                 
   for h=1:1:sz            
     MS(h,1)= abs(MC(op,1)-MC(h,1)); % take difference of present and previous value in column matrix MC
     op=op+1;               
   end
end

% *****    END of DFT Standard deviation      ****** %

% *****    threshold choosen     ******

th_h=max(MS);
th_l=th_h./8;
            
            
sz=size(MS,1);
q=1;
 for p=1:1:sz
    if(MS(p,1)>=th_l)                     % threshold window bet th_h & th_l    
        if(MS(p,1)<th_h)
        frNO(q,1)= p;                     % Returns frame number where change has occured
         q=q+1;
        end
     end
 end
  
frNO=frNO*step;



