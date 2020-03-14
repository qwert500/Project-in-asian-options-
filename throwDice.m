clc, clear all
N=100001; % number of tosses, e.g odd
rounds=1000;
firstTossIsHeads=0;
mostTossesIsHeads=0;

for j=1:rounds
  heads=0;
  
  %heads adds 1 to the sum, tails is 0
  
  for i=1:N
    x=randi([1 2],1,1);
    if x==2
      heads=heads+1;
    end
    if i==1
      if heads==1
        firstTossIsHeads=firstTossIsHeads+1;
      end
    end
    
  end
  
  if heads>N/2
    mostTossesIsHeads=mostTossesIsHeads+1;
  end
end

if firstTossIsHeads>rounds/2
  disp('most first tosses are heads')
end

if mostTossesIsHeads>rounds/2
  disp('most tosses are heads')
end


