%Script for generating plots based on color space regressions and make
%predictions based on those predictions

%clear all

%%% 0 is holdout, 1 is no holdout
mode = 1;

color_coords = readmatrix('../../data/UW58_colors.csv'); %%file with all color coordinates xyy Lab ch
color_coords = color_coords(:,4:6);
col = colorconvert( color_coords, 'Lab', 'D65');


if mode == 1

    for this_col = 1:58
    close all;
    clf;
    cors=zeros(2,8)
    predictions = zeros(8,58)
    C = cell(40,1); y = cell(40,1);  % 40 = upper bound on number of expts
    A = cell(30,1); L = cell(30,1);  % 30 = upper bound on number of models
    d = 1;
    r=0;
    c=0;
    l=1;

    %col = colorconvert( color_coords, 'xyY', 'D65', 'BCPbg' );


    figure(30);
    num_pc=8;

    filename = strcat('../../data/8_basis_vecs_scaled_color_',string(this_col),'.csv');
    basis_vecs = readmatrix(filename,  'Range', [2,2]);
    all_weights = zeros(7,8)

    for pc = 1:8

        figure(pc);
        L_weights=[]
        C_weights=[]
        r=0;
        ndatasets = 1;
        C = col;
        data = basis_vecs(:,pc)
        ndatasets = 1;

        y = data;
        expt_name = string(pc);
        On = ones(size(y,1),1);
        achrom = zeros(size(y,1),1);
        achrom(23:27) = 1;
        i=1;
        % regressors used (constant, L, 1st harmonic, 2nd harmonic, and Cab)
        C.L(this_col) = [];
        C.hab(this_col)=[];
        C.Cab(this_col)=[];
        A{i}=[ On C.L cosd(C.hab) sind(C.hab) cosd(2*C.hab) sind(2*C.hab) C.Cab];

        L{i}= 'LChab polar, 2nd'; 
        i=i+1;
        %% evaluate each model
        nmodels = i-1;
        [ncolors,nsubjects] = size(y);
        % compute weights for this set of regressors.
        for i = 1:nmodels
            weights = A{i}\y;
        end

        cor =  corrcoef(y,(A{1}*weights));

        %cors = [cors,cor(1,2)];
        %cors(r, c) = cor(1,2);
        L_weights = [L_weights,weights(2)];
        C_weights = [C_weights,weights(7)];
        all_weights(:,pc) = weights;

        %% convert weights into dominant hue angles

        % hue angle is given by cos and sin. Convert to angle + magnitude
        hangle = mod( atan2d( weights(4,:), weights(3,:) ), 360 )';
        rho = sqrt(sum(weights(3:4,:).^2,1))';

        % same thing for second harmonic. Take angle as well as angle+180.
        hangle2 = mod( [ atan2d( weights(6,:), weights(5,:) )'/2;
                         atan2d( weights(6,:), weights(5,:) )'/2 + 180], 360 );
        rho2 = sqrt(sum(weights(5:6,:).^2,1))'; rho2 = [rho2;rho2];

        % convert hue angles to RGB so we can color the plots.

        FIT = AllTheColors( [78+14*sind(hangle) 60*ones(nsubjects,1) hangle], 'Lchab', 'D65', 'D65' );
        [R,G,B] = Lab2RGB( FIT.L, FIT.a, FIT.b );
        cols = [R G B];

        FIT2 = AllTheColors( [78+14*sind(hangle2) 60*ones(2*nsubjects,1) hangle2], 'Lchab', 'D65', 'D65' );
        [R2,G2,B2] = Lab2RGB( FIT2.L, FIT2.a, FIT2.b );
        cols2 = [R2 G2 B2];

% 
%         f = subplot(1,1,1);
%         for i = 1:nsubjects
% 
%             polarplot( [hangle(i);hangle(i)]*pi/180, [0;rho(i)], 'Color', cols(i,:),...
%                 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 1.5, 'LineStyle', ':', 'MarkerFaceColor', cols(i,:) )
%             text(hangle(i)*pi/180,rho(i)+10, string(pc))
%             hold on;
% 
%         end
%         title([strcat('PC',expt_name) ', 1st harmonic'])
%         rlim([0 0.1])
% 
%         f = subplot(1,1,1);
%         f.RTick = [30 80];
%         f.ThetaTick = 0:45:360;
%         f.FontSize = 9;
% 
% 
%         for i = 1:2*nsubjects
% 
%             polarplot( [hangle2(i);hangle2(i)]*pi/180, [0;rho2(i)], 'Color', cols2(i,:),...
%                 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 1, 'MarkerFaceColor','w', 'LineStyle',':' )
%             text(hangle2(i)*pi/180,rho2(i)+10,string(pc))
%             hold on;          
%         end
%         title([strcat('PC',expt_name) ', 1st and 2nd Harmonics'])
%         rlim([0 0.1])
%         f.RTick = [30 80];
%         f.ThetaTick = 0:45:360;
%         f.FontSize = 9;
% 
% 
% 
% 
% 
%         ff = figure(pc);
%         ff.Position = [2 32 1000 700]; 
% 


    %     predictions(pc,:) = transpose(A{1}*weights);
    %     csvwrite('../../data/regression_predictions.csv',predictions)
    %      csvwrite('../../data/uw58_regressor_vals.csv',transpose(A{1}))
    %     

    %csvwrite('../../data/8_regression_weights.csv',all_weights)
    



    end
    filename = strcat('../../data/8_regression_weights_color_',string(this_col),'.csv');
    csvwrite(filename,all_weights);  

    end
    
    
    %--------------%
elseif mode == 0
    
    close all;
    clf;
    cors=zeros(2,8)
    predictions = zeros(8,58)
    C = cell(40,1); y = cell(40,1);  % 40 = upper bound on number of expts
    A = cell(30,1); L = cell(30,1);  % 30 = upper bound on number of models
    d = 1;
    r=0;
    c=0;
    l=1;


  
    num_pc=8;
    %figure(1)

    basis_vecs = readmatrix('../../data/8_basis_vecs_scaled.csv',  'Range', [2,2]);
    all_weights = zeros(7,8)
    
    
    
    for pc = 1:8
            L_weights =[]
            C_weights=[]

        figure(pc);
        
        r=0;
        ndatasets = 1;
        C = col;
        data = basis_vecs(:,pc)
        ndatasets = 1;

        y = data;
        expt_name = string(pc);
        On = ones(size(y,1),1);
        achrom = zeros(size(y,1),1);
        achrom(23:27) = 1;
        i=1;
        % regressors used (constant, L, 1st harmonic, 2nd harmonic, and Cab)
        A{i}=[ On C.L cosd(C.hab) sind(C.hab) cosd(2*C.hab) sind(2*C.hab) C.Cab];

        L{i}= 'LChab polar, 2nd'; 
        i=i+1;
        %% evaluate each model
        nmodels = i-1;
        [ncolors,nsubjects] = size(y);
        % compute weights for this set of regressors.
        for i = 1:nmodels
            weights = A{i}\y;
        end

        cor =  corrcoef(y,(A{1}*weights));

        %cors = [cors,cor(1,2)];
        %cors(r, c) = cor(1,2);
        L_weights = [L_weights,weights(2)];
        C_weights = [C_weights,weights(7)];
        all_weights(:,pc) = weights;

        %% convert weights into dominant hue angles

        % hue angle is given by cos and sin. Convert to angle + magnitude
        hangle = mod( atan2d( weights(4,:), weights(3,:) ), 360 )';
        rho = sqrt(sum(weights(3:4,:).^2,1))';

        % same thing for second harmonic. Take angle as well as angle+180.
        hangle2 = mod( [ atan2d( weights(6,:), weights(5,:) )'/2;
                         atan2d( weights(6,:), weights(5,:) )'/2 + 180], 360 );
        rho2 = sqrt(sum(weights(5:6,:).^2,1))'; rho2 = [rho2;rho2];

        % convert hue angles to RGB so we can color the plots.

        FIT = AllTheColors( [78+14*sind(hangle) 60*ones(nsubjects,1) hangle], 'Lchab', 'D65', 'D65' );
        [R,G,B] = Lab2RGB( FIT.L, FIT.a, FIT.b );
        cols = [R G B];

        FIT2 = AllTheColors( [78+14*sind(hangle2) 60*ones(2*nsubjects,1) hangle2], 'Lchab', 'D65', 'D65' );
        [R2,G2,B2] = Lab2RGB( FIT2.L, FIT2.a, FIT2.b );
        cols2 = [R2 G2 B2];


        f = subplot(1,1,1);
        for i = 1:nsubjects

            polarplot( [hangle(i);hangle(i)]*pi/180, [0;rho(i)], 'Color', cols(i,:),...
                'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 1.5, 'LineStyle', '-', 'MarkerFaceColor', cols(i,:) )
            text(hangle(i)*pi/180,rho(i)+10, string(pc))
            hold on;

        end
        title([strcat('PC',expt_name) ', 1st harmonic'])
        rlim([0 0.1])

        f = subplot(1,1,1);
        f.RTick = [30 80];
        f.ThetaTick = 0:45:360;
        f.FontSize = 15;


        for i = 1:2*nsubjects

            polarplot( [hangle2(i);hangle2(i)]*pi/180, [0;rho2(i)], 'Color', cols2(i,:),...
                'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 1, 'MarkerFaceColor','w', 'LineStyle',':' )
            text(hangle2(i)*pi/180,rho2(i)+10,string(pc))
            hold on;          
        end
        title([strcat('PC',expt_name) ', 1st and 2nd Harmonics'])
        rlim([0 0.15])
        f.RTick = [.05 .10];
        f.ThetaTick = 0:45:360;
        f.FontSize = 15;





        ff = figure(pc);
        ff.Position = [2 32 1000 700]; 



    %     predictions(pc,:) = transpose(A{1}*weights);
    %     csvwrite('../../data/regression_predictions.csv',predictions)
    %      csvwrite('../../data/uw58_regressor_vals.csv',transpose(A{1}))
    %     

    
    print(ff,strcat('dand_plots_VSS/dandelion_grid_',string(pc)),'-dpdf')

    
fig = figure(30);
f_ = subplot(2,num_pc,l);
%plot([-2 2],L_weights, 'Marker','.','MarkerSize',20);
bar(L_weights);
ylim([-0.005 0.005]);
title([strcat('PC',string(pc), '\_lightness'),''])
ylabel('Weight')
xlabel('')
xticks('')
hold on;




f_ = subplot(2,num_pc,num_pc+l);
%plot([-2 2],C_weights, 'Marker','.','MarkerSize',20);
bar(C_weights);
ylim([-0.005 0.005]);
title([strcat('PC',string(pc), '\_chroma'),''])
ylabel('Weight')
xlabel('')
xticks('')
hold on;

l=l+1;
    
predictions(pc,:) = transpose(A{1}*weights);


    end
    csvwrite('../../data/8_regression_weights.csv',all_weights)
    orient(fig,'landscape')
    set(fig,'PaperUnits','inches','PaperPosition',[0 0 15 6], 'PaperSize',[16 7]);
    print(fig,'dand_plots_VSS/lightness_chroma','-dpdf');
    csvwrite('../../data/regression_predictions.csv',predictions)
    

end



% 
% for pc = 1:8
% %clf;
% c=c+1;
% figure(pc); 
% L_weights=[]
% C_weights=[]
% r=0;
% for sd = [-2 2]
%     r=r+1;
%     filename = strcat('../../data/recon_',string(pc),'_',string(sd),'.csv');
%     data = readmatrix(filename,'NumHeaderLines',1);
%     data = transpose(data);
%     %data = rescale(data, -100,100);
%     
%     ndatasets = 1;
% 
%     C = col;
%     y = data;
%     expt_name = string(pc);
%     On = ones(size(y,1),1);
%     achrom = zeros(size(y,1),1);
%     achrom(23:27) = 1;
%     i=1;
%     % regressors used (constant, L, 1st harmonic, 2nd harmonic, and Cab)
%     A{i}=[ On C.L cosd(C.hab) sind(C.hab) cosd(2*C.hab) sind(2*C.hab) C.Cab];
%     
%     
%     L{i}= 'LChab polar, 2nd'; 
%     i=i+1;
% 
%     %% evaluate each model
%     nmodels = i-1;
%     [ncolors,nsubjects] = size(y);
% 
%     % compute weights for this set of regressors.
%     for i = 1:nmodels
%         weights = A{i}\y;
%     end
%     
%     cor =  corrcoef(y,(A{1}*weights));
%     disp(weights)
%     %cors = [cors,cor(1,2)];
%     cors(r, c) = cor(1,2);
%     L_weights = [L_weights,weights(2)];
%     C_weights = [C_weights,weights(7)];
% 
%     %% convert weights into dominant hue angles
% 
%     % hue angle is given by cos and sin. Convert to angle + magnitude
%     hangle = mod( atan2d( weights(4,:), weights(3,:) ), 360 )';
%     rho = sqrt(sum(weights(3:4,:).^2,1))';
% 
%     % same thing for second harmonic. Take angle as well as angle+180.
%     hangle2 = mod( [ atan2d( weights(6,:), weights(5,:) )'/2;
%                      atan2d( weights(6,:), weights(5,:) )'/2 + 180], 360 );
%     rho2 = sqrt(sum(weights(5:6,:).^2,1))'; rho2 = [rho2;rho2];
% 
%     % convert hue angles to RGB so we can color the plots.
% 
%     FIT = AllTheColors( [78+14*sind(hangle) 60*ones(nsubjects,1) hangle], 'Lchab', 'D65', 'D65' );
%     [R,G,B] = Lab2RGB( FIT.L, FIT.a, FIT.b );
%     cols = [R G B];
% 
%     FIT2 = AllTheColors( [78+14*sind(hangle2) 60*ones(2*nsubjects,1) hangle2], 'Lchab', 'D65', 'D65' );
%     [R2,G2,B2] = Lab2RGB( FIT2.L, FIT2.a, FIT2.b );
%     cols2 = [R2 G2 B2];
% 
%   
% 
%     %f = subplot(2,num_pc,pc);
%     f = subplot(1,1,1);
%     
%     for i = 1:nsubjects
%         if sd>0
%             polarplot( [hangle(i);hangle(i)]*pi/180, [0;rho(i)], 'Color', cols(i,:),...
%                 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 4, 'MarkerFaceColor', cols(i,:) )
%             text(hangle(i)*pi/180,rho(i)+10,string(sd))
%             hold on;
%         else
%             polarplot( [hangle(i);hangle(i)]*pi/180, [0;rho(i)], 'Color', cols(i,:),...
%                 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 1.5, 'LineStyle', ':', 'MarkerFaceColor', cols(i,:) )
%             text(hangle(i)*pi/180,rho(i)+10,string(sd))
%             hold on;
%         end
%     end
%     title([strcat('PC',expt_name) ', 1st harmonic'])
%     rlim([0 0.3])
%     
%     %f = subplot(2,num_pc,pc);
%     f = subplot(1,1,1);
% 
%     f.RTick = [30 80];
%     f.ThetaTick = 0:45:360;
%     f.FontSize = 9;
%     
% 
%     %f = subplot(2,num_pc,pc+num_pc);
%     
%     %f = subplot(2,1,2);
%     
%     for i = 1:2*nsubjects
%         if sd>0
%             polarplot( [hangle2(i);hangle2(i)]*pi/180, [0;rho2(i)], 'Color', cols2(i,:),...
%                 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 1, 'MarkerFaceColor','w' )
%             text(hangle2(i)*pi/180,rho2(i)+10,string(sd))
%             hold on;
%         else
%              polarplot( [hangle2(i);hangle2(i)]*pi/180, [0;rho2(i)], 'Color', cols2(i,:),...
%                 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 1, 'MarkerFaceColor','w', 'LineStyle',':' )
%             text(hangle2(i)*pi/180,rho2(i)+10,string(sd))
%             hold on;
%         end
%             
%     end
%     %title([strcat('PC',expt_name) ', 2nd harmonic'])
%     
%     title([strcat('PC',expt_name) ', 1st and 2nd Harmonics'])
%     rlim([0 0.3])
%     
%     %f = subplot(2,num_pc,pc+num_pc);
%     %f = subplot(2,1,2);
%     
%     f.RTick = [30 80];
%     f.ThetaTick = 0:45:360;
%     f.FontSize = 9;
% 
% 
%     %     figure(2)
%     %     subplot(1,ndatasets,d)
%     %     for i = 1:nsubjects
%     %         plot( hangle(i), hangle(i)-hangle2(i+nsubjects), 'Color', cols(i,:),...
%     %             'Marker', '.', 'MarkerSize', 20,'LineStyle', 'none' )
%     %         hold on;
%     %         plot( hangle(i), hangle(i)-hangle2(i), 'Color', cols(i,:),...
%     %             'Marker', '.', 'MarkerSize', 20,'LineStyle', 'none' )
%     %     end
%     %     ylim([-90,90])
% 
% 
% 
%     ff = figure(pc);
%     %ff = figure(1);
%     ff.Position = [2 32 1000 700];
%   
%    
% 
%     
% if sd == -2
%     predictions(pc,:) = transpose(A{1}*weights);
%     csvwrite('../../data/regression_predictions.csv',predictions)
% end
% 
% 
% end


