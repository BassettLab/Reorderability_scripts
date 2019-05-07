%%
close all
clear

filenames = {'NF_ER_p04_1218','NF_PropProb_d1_1218','NF_abssin_2pi_1218',...
    'NF_prefAttach_m4_1218','NF_RG_ep015_1218','NF_spatialGrowth_b1a4_1218'};
for ngraph = 1:6
    filename = filenames{ngraph};
    load(sprintf('%s_simratiocomms.mat',filename))

    figure('Color',[1 1 1],'Position',[0 0 1200 900])

    % Define colormap 
    cmap1 = (1/256)*[linspace(256,93,1000)' linspace(256,113,1000)' linspace(256,178,1000)'];

    subplot(4,1,1)

    imagesc(node_sim_strength);
    colormap(gca,cmap1)
    colorbar


    subplot(4,1,2)

    % Scale node_sim_strength to match cmap1
    nssn = -node_sim_strength;
    cmap2a = 93+ [(nssn-min(nssn))./(max(nssn) - min(nssn))].*(256-93);
    cmap2b = 113+ [(nssn-min(nssn))./(max(nssn) - min(nssn))].*(256-113);
    cmap2c = 178+ [(nssn-min(nssn))./(max(nssn) - min(nssn))].*(256-178);
    cmap2 = (1/256)*[cmap2a cmap2b cmap2c];


    plot([0 70], [1 1],'w')
    hold on
    counter = 1;
    for i = 1:70
        for j = (i+1):70
            if badj(i,j) == 1

                plotHalfCircle([mean([i j]),1],(abs(i-j)/2));
                counter = counter+1;
            end
        end
        % Plot node
        plot(i,1,'.','MarkerSize',16,'Color',cmap2(i,:))
    end

    hold off

    axis equal
    axis off





    subplot(4,1,3)

    % Generate colormap2
    cmap3 = rand(max(simRatio_comms),3);
    disp(max(simRatio_comms))
    plot([0 70], [1 1],'w')
    hold on
    counter = 1;
    for i = 1:70
        for j = (i+1):70
            if badj(i,j) == 1

                plotHalfCircleNeg([mean([i j]),1],(abs(i-j)/2));
                counter = counter+1;
            end
        end
        plot(i,1,'.','Color',cmap3(simRatio_comms(i),:),'MarkerSize',16)
    end
    hold off

    axis equal
    axis off


    subplot(4,1,4)
    image(simRatio_comms)
    title(sprintf('Modularity = %f, nComms = %i',Q,max(simRatio_comms)))
    colormap(gca,cmap3)
    colorbar
    print(sprintf('%s_circleplot_comms',filename),'-dpdf','-painters')

    disp('Done :)')
end