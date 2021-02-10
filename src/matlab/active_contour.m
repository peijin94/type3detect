function real_line = active_contour(t,f,sub_img,init_point_a,init_point_b,n_point)
    %init xy serise
    iter_len = 901;
    delta_go = 0.003;
    alpha=25;
    beta=75;
    [T,F] = meshgrid(t,f);
    rem_ends = ones(1,n_point);
    rem_ends(1)=-2.1;
    rem_ends(end)=-2.1;
    
    xx = linspace(init_point_a(1),init_point_b(1),n_point);
    yy = linspace(init_point_a(2),init_point_b(2),n_point)*70/400+10;
    
    verts = [xx;yy]';
    
    curv = line_curvature(verts);
    
    
    % process sub_img for better result
    
    sub_img_tmp = sub_img;
    
    h0=fspecial('motion',12,90);
    h1=fspecial('disk',3);
    img_o = imfilter(imfilter(sub_img_tmp,h0,'replicate'),h1,'replicate'); 
    sub_img=img_o;
    
    [px,py] = gradient((sub_img)/max(sub_img(:)));
    %quiver(t,f,px,py)
    interp2(T,F,px,xx,yy,'cubic');

    
    xx_new=xx;
    yy_new=yy;
    
    
    for num=1:iter_len
    
        xx_new = xx_new + delta_go*(alpha* interp2(T,F,px,xx_new,yy_new,'cubic') -rem_ends.* beta.*line_curvature([xx_new;yy_new]')'); 
        yy_new = yy_new + delta_go/12*(alpha* interp2(T,F,py,xx_new,yy_new,'cubic') -rem_ends.* beta.*line_curvature([xx_new;yy_new]')'); 
        
%         hp1=plot(xx_new,yy_new,'wo-','linewidth',3,'markersize',3);
%         hp2=plot(xx_new,yy_new,'ko-','linewidth',2,'markersize',1.5);
%         if mod(num-1,10)==0
%             set(hp1,'visible','on')
%             set(hp2,'visible','on')
%             ht=text(169,75,['iter : ',num2str(num-1)],'color','w','fontsize',19);
%              set(gca,'position',[0 0 1 1])
%              set(gcf,'position',[186   405   255   316])
%             xlim([155 178])
%             saveas(gcf,['img/act/act_',num2str((num-1)),'.eps'],'epsc')
%             
%             set(hp1,'visible','off')
%             set(hp2,'visible','off')
%         end
%         delete(hp1)
%         delete(hp2)
%         delete(ht)
    end
    hp1=plot(xx_new,yy_new,'wo-','linewidth',3,'markersize',3);hold on
    hp2=plot(xx_new,yy_new,'ko-','linewidth',2,'markersize',1.5);
    
            set(hp1,'visible','on')
            set(hp2,'visible','on')
    real_line.xx=xx_new;
    real_line.yy=yy_new;
end