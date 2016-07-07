function heatmap_single_bump_inR(R)
hw = R.grid.hw;
fw = 2*hw+1;
Mat = zeros(fw);
time = datestr(now,'yyyymmddTHHMMSS');
if size(R.grid.centre) == [0,0]
    imagesc(Mat)
else
    x_mean_chosen = R.grid.centre(1,:);
    y_mean_chosen = R.grid.centre(2,:);
    dist_std_chosen = R.grid.radius;
    t_mid = R.grid.t_mid_full;
    t_mid_chosen = R.grid.t_mid;
    j = 1;
    for t = 1:length(t_mid);
        if sum(t_mid_chosen == t_mid(t)) == 1
            x_tmp = round(x_mean_chosen(j)) + hw + 1;
            y_tmp = round(y_mean_chosen(j)) + hw + 1;
            Mat(x_tmp,y_tmp) = Mat(x_tmp,y_tmp) + 1/dist_std_chosen(j);
            j = j + 1;
        end
    end
    Mat = permute(Mat,[2 1]);
    Mat = flipud(Mat);
    imagesc(Mat)
end
colorbar
tit = sprintf('loop number = %04i time:%s',R.ExplVar.loop_num,time);
name = sprintf('%04i_heatmap_%s.pdf',R.ExplVar.loop_num,time);
title(tit)
set(gcf,'renderer','zbuffer');
saveas(gcf,name);
end