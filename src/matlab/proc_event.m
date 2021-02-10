function proc_event(t,f,img,img0,lines1,lines2,outdir,fname)
    sz_img = size(img);
    [lines1_neat,xypoint_1,xypoint_1h]=lines_neater(lines1);
    [~,xypoint_2,xypoint_2h]=lines_neater(lines2);
    
    index_lines = 1:length(lines1_neat);
    index_pool  = index_lines;
    
    idx_num=0;
    
    % split the events in bins
    while ~isempty(index_pool)
        idx_num=idx_num+1;
        index_select=index_pool(1);
        index_in_pool = ...
            (abs(xypoint_1h(index_pool,1)...
            -xypoint_1h(index_select,1))<=9);

        event_set(idx_num).indexs = index_pool(index_in_pool);
        index_pool(index_in_pool)=[];
        event_set(idx_num).range_x = [min(xypoint_1h(event_set(idx_num).indexs,1))...
            ,max(xypoint_1h(event_set(idx_num).indexs,1))];
        event_set(idx_num).bad_e_flag=0;
        
        if event_set(idx_num).range_x(1)>268 || event_set(idx_num).range_x(2)<32
             event_set(idx_num).bad_e_flag=1;
        end
    end
    
    for num=1:length(event_set)
        event_cur = event_set(num);
        % if not a bad event
        if ~event_cur.bad_e_flag
            focus_range = [max(event_cur.range_x(1)-20,1),min(event_cur.range_x(2)+20,sz_img(2))];        
            sub_img = img(:,focus_range(1):focus_range(2));
            sub_img0 = img0(:,focus_range(1):focus_range(2)); 
            
            for num_i=1:length(event_cur.indexs)
                %line_cur = lines1(event_cur.indexs(num));
            %    plot([line_cur.point1(1),line_cur.point2(1)]...
            %        ,[line_cur.point1(2),line_cur.point2(2)]*70/400+10,'w','linewidth',2);
            end

            esti_low = [mean(mean(xypoint_1(event_cur.indexs,1))),...
                mean(mean(xypoint_1(event_cur.indexs,2)))];
            fake_high = [mean(mean(xypoint_1h(event_cur.indexs,1))),...
                mean(mean(xypoint_1h(event_cur.indexs,2)))];

            drift_range = mean(mean(xypoint_1h(event_cur.indexs,2)))...
                -mean(mean(xypoint_1(event_cur.indexs,2)));
            % find lines to estimate the slope
            estimate_idx = find(xypoint_2h(:,1)>focus_range(1)-10 & xypoint_2h(:,1)<10+focus_range(2));
            if ~isempty(estimate_idx)
                esti_slope = mean(xypoint_2h(estimate_idx,2)-xypoint_2(estimate_idx,2))/ ...
                    mean(xypoint_2h(estimate_idx,1)-xypoint_2(estimate_idx,1));

                esti_hx = esti_low(1)+(drift_range)/esti_slope;
                esti_hy = mean(mean(xypoint_1h(event_cur.indexs,2)));


                h=figure();

                imagesc(t(focus_range(1):focus_range(2)),f,(sub_img));
                
                set(gca,'clim',[prctile(sub_img(:),1.5),prctile(sub_img(:),98.5)])
                colormap(jet)
                hold on
                xlabel('Time (s)')
                ylabel('Frequency (MHz)')
                esti_high=[esti_hx,esti_hy];
                line_esti_high = [esti_high(1)-(esti_high(1)-fake_high(1)),esti_high(2)];        
                line_esti_low = [esti_low(1)-(esti_high(1)-fake_high(1)),esti_low(2)];
                plot([esti_high(1),esti_low(1)]-(esti_high(1)-fake_high(1)),[esti_high(2),esti_low(2)]*70/400+10,'w','linewidth',2)
                real_line = active_contour(t(focus_range(1):focus_range(2)),f,sub_img,line_esti_low,line_esti_high,16);

                [~,dtstr,idtime] = fname2datetime(fname);
                date_str = dtstr(1:10);
                time_str = dtstr(12:end);
                if ~exist([outdir,'/',date_str(1:7)],'dir')
                    mkdir([outdir,'/',date_str(1:7)]);
                end
                print(h,[outdir,'/',date_str(1:7),'/[',date_str(9:10),']',replace(time_str,':','_'),'_event','_id_',idtime,'[',num2str(num),'].jpg'],'-djpeg','-r300')
                %close(h)

            
                
                h_act=figure();
                subplot(121)
                imagesc(t(focus_range(1):focus_range(2)),f,(sub_img0));
                colormap(jet)
                hold on
                plot([esti_high(1),esti_low(1)]-(esti_high(1)-fake_high(1)),[esti_high(2),esti_low(2)]*70/400+10,'w','linewidth',2)
                xx_new=real_line.xx;
                yy_new=real_line.yy;
                
                hp1=plot(xx_new,yy_new,'wo-','linewidth',3,'markersize',3);hold on
                hp2=plot(xx_new,yy_new,'ko-','linewidth',2,'markersize',1.5);
            
                t_base=linspace(xx_new(1),xx_new(end),200);
                f_base=interp1(xx_new,yy_new,t_base,'spline');
                
                dfdt=diff(f_base)/(t_base(2)-t_base(1));
                
                xlabel('Time (s)')
                ylabel('Frequency (MHz)')
                
                subplot(222)
                
                %plot(f_base(1:(end-1)),dfdt,'k-s')
                plot((yy_new(1:end-1)+yy_new(2:end))/2,diff(yy_new)./(xx_new(2:end)-xx_new(1:end-1)),'k-s')
                ylabel('Drift velocity (MHz/s)')
                xlim([f_base(1),f_base(end)])
                set(gca,'position',[0.5703    0.520    0.3347    0.3900])
                ylim([-13 0])
                legend('Drift velocity')
                set(gca,'xtick',[])
                subplot(224)
                
                [T,F] = meshgrid(t,f);
                line0 = interp2(T,F,img0,t_base,f_base,'cubic');
                line  = interp2(T,F,img,t_base,f_base,'cubic');
                
                
                
                yyaxis left
                plot(f_base,line0)
                ylabel('Original')
                ylim([72,135])
                yyaxis right
                plot(f_base,line,'--')
                ylim([0,60])
                ylabel('Preprocessed')
                
                xlim([f_base(1),f_base(end)])
                set(gca,'position',[0.5703    0.130    0.3347    0.3900])
                set(gcf,'position',[413.0000  439.4000  678.4000  403.2000])
                xlabel('Frequency (MHz)')
                legend('Original','Preprocessed')
                saveas(h_act,[outdir,'/',date_str(1:7),'/[',date_str(9:10),']',replace(time_str,':','_'),'_more','_id_',idtime,'[',num2str(num),'].eps'],'epsc')
                
%                print(h_act,[outdir,'/',date_str(1:7),'/[',date_str(9:10),']',replace(time_str,':','_'),'_more','_id_',idtime,'[',num2str(num),'].jpg'],'jpeg','-r300')
                t_get=t(focus_range(1):focus_range(2));
                    f_get=f;
                    img_get=sub_img;
                    save([outdir,'/meta/[',date_str(1:4),date_str(6:7),date_str(9:10),']',...
                        replace(time_str,':','_'),'_event','_id_',idtime,'[',num2str(num),'].mat'],...
                        't_get','f_get','img_get','real_line','fname')

                
            end
        end
    end 
    
end