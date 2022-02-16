% runMe

% define parameters
fill_gap_len = 5;
min_len_along = 56;
thresh_cap = 0.86;
hough_peak1_thresh=270;
hough_peak2_thresh=80;
nhoodsize=[9,5];

% normal 91-120
hough_angle_lower = 91;
hough_angle_upper = 101;
hough_angle_res=0.2;
t=1:300;
f=(1:400)/400*70+10;
T_range=0.9:0.1:11.5;

%
% % sample file used for debug 
% % fname0 = '../data/20130103_l_38645.fits';
% % fname1 = '../data/20130113_l_47285.fits';
% % fname2 = '../data/20130326_l_45905.fits';
% % fname3 = '../data/20131019_l_41345.fits';
% 

outdir='outdir_3/';
names = ls('I:\database\nancay\process\src\processed_tmp\*l*.fits');
tic;
for i=1:1000%length(names(:,1))
    name = ['I:\database\nancay\process\src\processed_tmp\',names(i,:)];
    name_l = ['I:\database\nancay\process\src\processed_tmp\',names(i,:)];
    name_r = strrep(['I:\database\nancay\process\src\processed_tmp\',names(i,:)],'r','l');
    %draw_normal(t,f,name_l,name_r,[outdir,'/all']);
    
    img0 = ((fitsread(name)'));
    img = constback_sub(img0);
    max_local = get_local_max_map(img);
    thresh_img = find(diff((hist(img(:),255)*triu(ones(255),0)>thresh_cap*length(img(:)))));
    BW = (img-thresh_img)>0;
    se = strel('ball',3,3);
    e_BW = edge(imerode(double(BW),se),'canny');
    [H,T,R] = hough(BW,'RhoResolution',2,'Theta',T_range);  % perform hough transform
    
    [H_edge,T_edge,R_edge] = hough(boolean(e_BW),'RhoResolution',2,'Theta',-0.4:0.1:0.4);
    run_flag=1;
    P_test  = houghpeaks(H_edge,6,'threshold',235);
    if ~isempty(P_test) % if event find!
        lines_t = houghlines(boolean(e_BW),T_edge,R_edge,P_test,'FillGap',fill_gap_len+3,'MinLength',min_len_along-5);
        acu_length=0;
        for num_t=1:length(lines_t)
            %length_line(lines_t(num_t).point1,lines_t(num_t).point2)
            if ... %length_line(lines_t(num_t).point1,lines_t(num_t).point2)>min_len_along && ...
                abs(lines_t(num_t).theta)<0.5%lines_t(num_t).point1(1)==lines_t(num_t).point2(1)
                acu_length=acu_length+length_line(lines_t(num_t).point1,lines_t(num_t).point2);
            end
        end
        if acu_length>350
            run_flag=0;
        end
    end
    H_ori=H;% save original H
%     H(:,1:hough_angle_lower)=0;
%     H(:,hough_angle_upper:180)=0;  % shield useless lines
%     
%     figure('visible','on')
%     imagesc(H)
    P  = houghpeaks(H,6,'Nhoodsize',nhoodsize,'threshold',hough_peak1_thresh);
    if ~isempty(P) && run_flag % if event find!
        % draw noise event
        % draw_noise_event(t,f,img,img0,BW,max_local,[outdir,'/','baddata/'],name_l);
        % be quiet or run: disp(['Something at :',name])
        
        lines = houghlines(BW,T,R,P,'FillGap',fill_gap_len,'MinLength',min_len_along);

        BW2 = max_local;

        [H2,T2,R2] = hough(BW2,'RhoResolution',2,'Theta',T_range);
        H2_ori=H2;
        %H2(:,1:hough_angle_lower)=0;
        %H2(:,hough_angle_upper:180)=0;

        P2  = houghpeaks(H2,6,'Nhoodsize',nhoodsize,'threshold',hough_peak2_thresh);
        if ~isempty(P2)&& ~isempty(lines)
            lines2 = houghlines(BW2,T2,R2,P2,'FillGap',ceil(fill_gap_len*1.5),'MinLength',min_len_along);
            if ~isempty(lines2)
            [lines_neat1,xypoint_a1,xypoint_a2]=lines_neater(lines);
            [lines_neat2,xypoint_b1,xypoint_b2]=lines_neater(lines2);
            
            ev_time = mean(xypoint_a2(:,1));
            
            dydx = (xypoint_b2-xypoint_b1);
            % [[x1 y1];[x2 y2];......[xn yn]]
            slope = (mean((dydx(:,2)/400*70)./mean(dydx(:,1))));
            
            [dt,dtstr,idtime]=fname2datetime(name);
            strtime = datestr(datenum(dt)+ev_time);
            strprint=['Event time : ', strtime,'      Frequency shift velo:',num2str(slope),' (MHz/s)'];
            disp(['[info] Type III event found !!!!  ',strprint]);
            draw_real_event(t,f,img,BW,max_local,H_ori,T,R,P,H2_ori,T2,R2,P2,lines,lines2,[outdir,'/events/'],name,strprint);
            
            proc_event(t,f,img,img0,lines_neat1,lines_neat2,[outdir,'/active_con'],name);
            end
        else
            % be quiet or disp('[bad data] ignore')
        end
        
    end
    close all
end
toc