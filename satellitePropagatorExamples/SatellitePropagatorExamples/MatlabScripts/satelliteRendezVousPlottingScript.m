close all
clear

PoweredChaserPhaseHistory = importdata('../chaserPropagationHistory_1.dat' );
TargetHistory = importdata('../targetPropagationHistory_1.dat' );

UnpoweredChaserPhaseHistory = importdata('../chaserPropagationHistory_2.dat' );
TargetHistory = [ TargetHistory; importdata('../targetPropagationHistory_2.dat' ) ];

j = 1;

%Translate to a planar view on the target orbital plane
for i = 1:length(TargetHistory(:,1))

   H = cross( TargetHistory(i, 2:4)', TargetHistory(i, 5:7)' );
   h(i,:) = (H/norm(H))';
   
   theta(i) = acos(dot(h(i,:)', [0;0;1]));
   
   E = cross( h(i,:)', [0;0;1]);
   
   e(i,:) = ( E/norm(E) )';
   
   quat = [ cos(theta(i)/2), e(i,1)*sin(theta(i)/2), e(i,2)*sin(theta(i)/2), e(i,3)*sin(theta(i)/2)];
   
   M_r = quat2rotm( quat );
   
   inPlaneTargetHistory(i,1:3) = (M_r*(TargetHistory(i, 2:4)'))';
   inPlaneTargetHistory(i,4:6) = (M_r*(TargetHistory(i, 5:7)'))';
   
   if i <= length(PoweredChaserPhaseHistory(:,1))
       
       inPlaneChaserPowered(i,1:3) = (M_r*(PoweredChaserPhaseHistory(i, 2:4)'))';
       inPlaneChaserPowered(i,4:6) = (M_r*(PoweredChaserPhaseHistory(i, 5:7)'))';
   
   else
       inPlaneChaserUnp(j,1:3) = (M_r*(UnpoweredChaserPhaseHistory(j, 2:4)'))';
       inPlaneChaserUnp(j,4:6) = (M_r*(UnpoweredChaserPhaseHistory(j, 5:7)'))';
       j = j+1;
   end
    
end


inPlaneTargetHistory(:,3) = zeros(length(inPlaneTargetHistory(:,3)), 1);

j=1;


%Translate to target-frame
for i = 1:length(inPlaneTargetHistory(:,1))
    
    angle = acos( dot( (inPlaneTargetHistory(i, 1:3 )'/norm(inPlaneTargetHistory(i, 1:3 ))), [1;0;0] ) );
    
    R = rotz( angle );
    
    if i <= length(PoweredChaserPhaseHistory(:,1))
        
        deltar = inPlaneChaserPowered(i,1:3) - inPlaneTargetHistory(i, 1:3);
        targetFramePoweredChaser(i,1:3) = (R*deltar')';       
   
    else
        
        deltar = inPlaneChaserUnp(j,1:3) - inPlaneTargetHistory(i, 1:3);
        targetFrameUnpoweredChaser(j,1:3) = (R*deltar')';   
        j = j+1;
    
    end
    
end

figure

plot( targetFramePoweredChaser(:,1), targetFramePoweredChaser(:,2), 'LineWidth', 4 );
hold on

plot( targetFrameUnpoweredChaser(:,1), targetFrameUnpoweredChaser(:,2) );

plot( 0, 0, '+' )

xlabel ('x (Radial) [m]')
ylabel ('y [m]')


axis equal

xlim( [-6e4, 14e4] );

grid on

legend( 'Chaser - Powered flight', 'Chaser - Coasting', 'Target' )


figure

plot( targetFramePoweredChaser(:,1), targetFramePoweredChaser(:,3), 'LineWidth', 4 );
hold on

plot( targetFrameUnpoweredChaser(:,1), targetFrameUnpoweredChaser(:,3) );

plot( 0, 0, '-+' )

axis equal

xlim( [-6e4, 14e4])

ylim( [-1e4, 1e4] )

xlabel( 'x (Radial diretion) [m]' )

ylabel( 'z [m]' ) 

grid on








