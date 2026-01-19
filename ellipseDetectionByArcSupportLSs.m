function [ellipses, L, posi] = ellipseDetectionByArcSupportLSs(I, Tac, Tr, specified_polarity)
%input:
% I: input image
% Tac: elliptic angular coverage (completeness degree)
% Tni: ratio of support inliers on an ellipse
%output:
% ellipses: N by 5. (center_x, center_y, a, b, phi)
% reference:
% 1��von Gioi R Grompone, Jeremie Jakubowicz, Jean-
% Michel Morel, and Gregory Randall, ��Lsd: a fast line
% segment detector with a false detection control.,�� IEEE
% transactions on pattern analysis and machine intelligence,
% vol. 32, no. 4, pp. 722�C732, 2010.
    angleCoverage = Tac;%default 165��
    Tmin = Tr;%default 0.6 
    unit_dis_tolerance = 2; %max([2, 0.005 * min([size(I, 1), size(I, 2)])]);%�ڵ��������̲�С��max(2,0.5%*minsize)
    normal_tolerance = pi/9; %�������̽Ƕ�20��= pi/9
    t0 = clock;
    if(size(I,3)>1)
        I = rgb2gray(I);
        [candidates, edge, normals, lsimg] = generateEllipseCandidates(I, 2, specified_polarity);%1,sobel; 2,canny
    else
        [candidates, edge, normals, lsimg] = generateEllipseCandidates(I, 2, specified_polarity);%1,sobel; 2,canny
    end
%    figure; imshow(edge);
%     return;
%    subplot(1,2,1);imshow(edge);%show edge image
%    subplot(1,2,2);imshow(lsimg);%show LS image
    t1 = clock;
    disp(['the time of generating ellipse candidates:',num2str(etime(t1,t0))]);
    candidates = candidates';%ellipse candidates matrix Transposition
    if(candidates(1) == 0)%��ʾû���ҵ���ѡԲ
        candidates =  zeros(0, 5);
    end
    posi = candidates;
    normals    = normals';%norams matrix transposition
    [y, x]=find(edge);%�ҵ���0Ԫ�ص���(y)����(x)������
%     ellipses = [];L=[];
%     return;
    [mylabels,labels, ellipses] = ellipseDetection(candidates ,[x, y], normals, unit_dis_tolerance, normal_tolerance, Tmin, angleCoverage, I);%���ĸ����� 0.5% 20�� 0.6 180�� 
    disp('-----------------------------------------------------------');
    disp(['running time:',num2str(etime(clock,t0)),'s']);
%     labels
%     size(labels)
%     size(y)
    warning('on', 'all');
     L = zeros(size(I, 1), size(I, 2));%����������ͼ��Iһ����С��0����L
     L(sub2ind(size(L), y, x)) = mylabels;%labels,���ȵ���edge_pixel_n x 1,�����i����Ե������ʶ���˵�j��Բ������б��Ϊj,����Ϊ0����С edge_pixel_n x 1;����ת���浽ͼ���У���ͼ���б��
%     figure;imshow(L==2);%LLL
%     imwrite((L==2),'D:\Graduate Design\��ͼ\edge_result.jpg');
end
%% ================================================================================================================================
%����1
%����
%candidates: ncandidates x 5
%points:     ��Ե���ص������(x,y),nx2,nΪ�ܹ��ı�Ե����
%lineLabels: ����Ӧ������(xi,yi)��ǣ���Ӧ������Ӧ���߶Σ�nx1,δ�����Ϊ0
%lines:      �߶β�����-B,A,xmid,ymid������(xmid,ymid)��Ӧ��Ӧ���߶��е㣬mx4��mΪ�ܹ�m���߶�
%���
%labels��    ���ȵ���n x 1,�����i����Ե������ʶ���˵�j��Բ������б��Ϊj,����Ϊ0����С n x 1
%C��   ʶ������ĶԳ����ģ������ᣬ�̰��ᣬ����ǣ�ÿһ�и�ʽ��(x,y,a,b,phi)
function [mylabels,labels, ellipses] = ellipseDetection(candidates, points, normals, distance_tolerance, normal_tolerance, Tmin, angleCoverage, E)
    labels = zeros(size(points, 1), 1);
    mylabels = zeros(size(points, 1), 1);%����
    ellipses = zeros(0, 5);
  
    %% ���������Ժܴ�ĺ�ѡ��Բ�������㼫���ϸ�Ҫ��ֱ�Ӽ�������SE(salient ellipses)��ͬʱ��~SE����goodness����pseudo order
    goodness = zeros(size(candidates, 1), 1);%��ʼ��ʱΪ0������⵽������Բʱֱ����ȡ����Ӧλ�õ�goodness(i) = -1��ǡ�
    for i = 1 : size(candidates,1)
        %ellipse circumference is approximate pi * (1.5*sum(ellipseAxes)-sqrt(ellipseAxes(1)*ellipseAxes(2))
        ellipseCenter = candidates(i, 1 : 2);
        ellipseAxes   = candidates(i, 3:4);
        tbins = min( [ 180, floor( pi * (1.5*sum(ellipseAxes)-sqrt(ellipseAxes(1)*ellipseAxes(2)) ) * Tmin ) ] );%ѡ����
        %ellipse_normals = computePointAngle(candidates(i,:),points);
        %inliers = find( labels == 0 & dRosin_square(candidates(i,:),points) <= 1 );  % +-1����������֧���ڵ�
        %���ټ��㣬ֻ������Բ��Ӿ����ڵı�Ե��(��Բ�еĳ���a>b),s_dx�洢�������points������
        s_dx = find( points(:,1) >= (ellipseCenter(1)-ellipseAxes(1)-1) & points(:,1) <= (ellipseCenter(1)+ellipseAxes(1)+1) & points(:,2) >= (ellipseCenter(2)-ellipseAxes(1)-1) & points(:,2) <= (ellipseCenter(2)+ellipseAxes(1)+1));
        inliers = s_dx(dRosin_square(candidates(i,:),points(s_dx,:)) <= 1);
        ellipse_normals = computePointAngle(candidates(i,:),points(inliers,:));
        p_dot_temp = dot(normals(inliers,:), ellipse_normals, 2); %���ٺ�ellipse_normals(inliers,:)��Ϊ���ٺ�ellipse_normals
        p_cnt = sum(p_dot_temp>0);%����֮�٣���һ�μ���ͳ�ƣ�����ΪC����ʱע�����Բʱ�ڵ㼫�Ե�ѡȡ����
        if(p_cnt > size(inliers,1)*0.5)
            %��������,Ҳ�����ں����    
            %ellipse_polarity = -1;
            inliers = inliers(p_dot_temp>0 & p_dot_temp >= 0.923879532511287 );%cos(pi/8) = 0.923879532511287, �н�С��22.5��  
        else
            %������ͬ,Ҳ�����ڰ���� 
            %ellipse_polarity = 1;
            inliers = inliers(p_dot_temp<0 & (-p_dot_temp) >= 0.923879532511287 );
        end
        inliers = inliers(takeInliers(points(inliers, :), ellipseCenter, tbins)); 
        support_inliers_ratio = length(inliers)/floor( pi * (1.5*sum(ellipseAxes)-sqrt(ellipseAxes(1)*ellipseAxes(2)) ));
        completeness_ratio = calcuCompleteness(points(inliers,:),ellipseCenter,tbins)/360;
        goodness(i) = sqrt(support_inliers_ratio*completeness_ratio); %goodness = sqrt(r_i * r_c)
        %{
        if( support_inliers_ratio >= Tmin && completeness_ratio >= 0.75 ) %300/360 = 0.833333333333333 and ratio great than Tmin
            goodness(i) = -1;
            if (size(ellipses, 1) > 0)
                s_flag = false;
                for j = 1 : size(ellipses, 1)
                    %��ʶ�������Բ���ܹ���֮ǰʶ�������Բ�ظ���pi*0.1 = 0.314159265358979
                    if (sqrt((ellipses(j, 1) - candidates(i, 1)) .^ 2 + (ellipses(j, 2) - candidates(i, 2)) .^ 2) <= distance_tolerance ...
                            && sqrt((ellipses(j, 3) - candidates(i, 3)) .^ 2 + (ellipses(j, 4) - candidates(i, 4)) .^ 2 ) <= distance_tolerance ...
                            && abs( ellipses(j, 5) - candidates(i, 5) ) <= 0.314159265358979) %pi/10 = 18��
                        s_flag = true;
                        labels(inliers) = j;%����ظ��ˣ��ͰѸñ�ǩת�Ƶ�֮ǰ��Բ���档
                        break;%������ѭ����������һ����ѭ��
                    end
                end
                if (~s_flag)%������ظ�������뵽ʶ���Բ(circles)��
                    labels(inliers) = size(ellipses, 1) + 1;
                    ellipses = [ellipses; candidates(i, :)];
                    %drawEllipses(candidates(i, :)',E);
                end
            else
                labels(inliers) = size(ellipses, 1) + 1;
                ellipses = [ellipses; candidates(i, :)];%���
                %drawEllipses(candidates(i, :)',E);
            end
        else
            goodness(i) = sqrt(support_inliers_ratio*completeness_ratio); %goodness = sqrt(r_i * r_c)
        end
        %}
    end
    %drawEllipses(ellipses',E);ellipses
    [goodness_descending, goodness_index] = sort(goodness,1,'descend');%here we can use pseudo order to speed up 
    candidates = candidates(goodness_index(goodness_descending>0),:);

    %%
%    t1 = clock;
    angles = [300; 210; 150; 90];%�ǶȴӴ�С��֤��������
    angles(angles < angleCoverage) = [];%ֻ��������angleCoverage�Ĳ���
    if (isempty(angles) || angles(end) ~= angleCoverage)%���angelsΪ���ˣ�����angles��С��~=angleCoverage�����angleCoverage�������
        angles = [angles; angleCoverage];
    end
%    disp('��ʼ��һ��Բ�������ǶȽ�����֤����ʼangleLoop����ÿ��ѭ��������һ��angleCoverage���Ժ�ѡԲ������֤������Բ���ϵ��ڵ����ͨ�Է��������������������ȷ������Ӷ��ҵ���ЧԲ��ͬʱ�޳���ЧԲ');
    for angleLoop = 1 : length(angles)
        idx = find(labels == 0);%labels��СΪ��Ե��������edge_nx1����ʼ��ʱlabelsȫΪ0���ҵ�labels�е���0������
        if (length(idx) < 2 * pi * (6 * distance_tolerance) * Tmin)%��idx����С��һ��ֵʱ
            break;
        end
        [L2, L, C, validCandidates] = subEllipseDetection( candidates, points(idx, :), normals(idx, :), distance_tolerance, normal_tolerance, Tmin, angles(angleLoop), E, angleLoop);
        candidates = candidates(validCandidates, :);%����logical����validCandidates�����޳�����������Բ��ʣ�µ�Բ����������һ��angleloop��֤
      % size(candidates)
      % disp(angleLoop)
        if (size(C, 1) > 0)
            for i = 1 : size(C, 1)
                flag = false;
                for j = 1 : size(ellipses, 1)
                    %��ʶ�������Բ���ܹ���֮ǰʶ�������Բ�ظ���pi*0.1 = 0.314159265358979
                    if (sqrt((C(i, 1) - ellipses(j, 1)) .^ 2 + (C(i, 2) - ellipses(j, 2)) .^ 2) <= distance_tolerance ...
                        && sqrt((C(i, 3) - ellipses(j, 3)) .^ 2 + (C(i, 4) - ellipses(j, 4)) .^ 2) <= distance_tolerance ...
                        && abs(C(i, 5) - ellipses(j, 5)) <= 0.314159265358979) %pi/10 = 18��
                        flag = true;
                        labels(idx(L == i)) = j;%����ظ��ˣ��ͰѸñ�ǩת�Ƶ�֮ǰ��Բ���档ע��ע�⣺idx�����������label�Ǳ�Ǹñ�Ե�������˵�j��Բ�ϣ�idx��labels����һά������(n x 1)��labels���Ե��points(n x 2)����һ��,��idx��һ��
                        %==================================================
                        mylabels(idx(L2 == i)) = j;
                        %==================================================
                        break;%������ѭ����������һ����ѭ��
                    end 
                end
                if (~flag)%������ظ�������뵽ʶ���Բ(circles)��
                    labels(idx(L == i)) = size(ellipses, 1) + 1;
                    %=================================================================
                    %%��ʾ��ϳ�Բʱ���õ��ڵ�  my code 
                    mylabels(idx(L2 == i)) = size(ellipses, 1) + 1;%����
                    %=================================================================
                    ellipses = [ellipses; C(i, :)];
                end
            end
        end
    end
%    t2 = clock;
%    disp(['�������֤ʱ�䣺',num2str(etime(t2,t1))]);
end


%% ================================================================================================================================
%����2
%����
%list��      �����ѡ��Բ�ĺͰ뾶��ϣ�(x,y,a,b,r)����С candidate_n x 5.
%points:     ��Ե���ص������(x,y),nx2,nΪ�ܹ��ı�Ե����
%normals:    ÿһ����Ե���Ӧ���ݶ�������normals��СΪnx2����ʽΪ(xi,yi)
 
%���
%labels��    �����i����Ե�����ڼ�⵽�˵�j��Բ����labels��i�и�ֵΪj������Ϊ0.������pointsһ�£�n x 1
%circles:    �˴μ�⵽��Բ,(x,y,z),����⵽detectnum�������СΪdetectnum x 3
%validCandidates: list�ĺ�ѡԲ�У������i��Բ����⵽�˻��߲�����Բ����(Բ�����ڵ���������)�����i��λ��Ϊfalse(��ʼ��ʱΪtrue)����������һ��angleloop�ִ���֤ʱ�����޳���������Ҫ�ظ���֤��
%                 validCandidates�Ĵ�СΪ candidate_n x 1.
function [mylabels,labels, ellipses, validCandidates] = subEllipseDetection( list, points, normals, distance_tolerance, normal_tolerance, Tmin, angleCoverage,E,angleLoop)
    labels = zeros(size(points, 1), 1);%��Ե���ص��������n,n x 1
    mylabels = zeros(size(points, 1), 1);%����
    ellipses = zeros(0, 5);
    ellipse_polarity = 0; %��Բ����
    max_dis = max(points) - min(points);
    maxSemiMajor = max(max_dis);%���Ŀ��ܰ뾶(�˴��ɸ�Ϊ/2)
    maxSemiMinor = min(max_dis);
    distance_tolerance_square = distance_tolerance*distance_tolerance;
    validCandidates = true(size(list, 1), 1);%logical��������С candidate_n x 1
    convergence = list;%��ѡ��Բ����
    for i = 1 : size(list, 1)
        ellipseCenter = list(i, 1 : 2);
        ellipseAxes = list(i, 3:4);
        ellipsePhi  = list(i,5);
        %ellipse circumference is approximate pi * (1.5*sum(ellipseAxes)-sqrt(ellipseAxes(1)*ellipseAxes(2))
        tbins = min( [ 180, floor( pi * (1.5*sum(ellipseAxes)-sqrt(ellipseAxes(1)*ellipseAxes(2)) ) * Tmin ) ] );%ѡ����
       
        %�ҵ����Բlist(i,:)��Բ���ϵ��ڵ�,find�����жϵĽ��Ϊlogical���������ڵ����Ӧ��points�е��ж�Ӧλ��Ϊ1������Ϊ0����СΪn x 1
        %ͨ��find�ҳ���0Ԫ�ص�������inliers���Ǳ�������ӦΪ�ڵ���е�ֵ������ inlier_n x 1
       
        %���д��������⣬δ���м��Է����������ڽ���2 * distance_tolerance�ļ����෴�Ĵ����ڵ㱻���룬�Ӷ���ϴ���.
%           inliers = find(labels == 0 & abs(sqrt((points(:, 1) - circleCenter(1)) .^ 2 + (points(:, 2) - circleCenter(2)) .^ 2) - circleRadius) <= 2 * distance_tolerance & radtodeg(real(acos(abs(dot(normals, circle_normals, 2))))) <= normal_tolerance);
        %===============================================================================================================================================================
        %���д����Ϊ���´���,my code
        
%         size(labels)
%         size(dRosin_square(list(i,:),points) )
%         ppp = (dRosin_square(list(i,:),points) <= 4*distance_tolerance_square)
%           if( i == 11 && angleCoverage == 165)
%          drawEllipses(list(i,:)',E);
%           end
        %�˴��ǳ���Ҫ������Ӧ��Ҫ��֮ǰ�ĸ���С�ķ�ΧѰ���ڵ�
        %ellipse_normals = computePointAngle(list(i,:),points);
        %inliers = find(labels == 0 & (dRosin_square(list(i,:),points) <= distance_tolerance_square) );  % 2.25 * distance_tolerance_square , 4*distance_tolerance_square.
        %���ټ��㣬ֻ������Բ��Ӿ����ڵı�Ե��(��Բ�еĳ���a>b),i_dx�洢���������points�е�����.
        i_dx = find( points(:,1) >= (ellipseCenter(1)-ellipseAxes(1)-distance_tolerance) & points(:,1) <= (ellipseCenter(1)+ellipseAxes(1)+distance_tolerance) & points(:,2) >= (ellipseCenter(2)-ellipseAxes(1)-distance_tolerance) & points(:,2) <= (ellipseCenter(2)+ellipseAxes(1)+distance_tolerance));
        inliers = i_dx(labels(i_dx) == 0 & (dRosin_square(list(i,:),points(i_dx,:)) <= distance_tolerance_square) );
        ellipse_normals = computePointAngle(list(i,:),points(inliers,:));%ellipse_normals������inliers����һ��
        
%         if( i == 11 && angleCoverage == 165)
%         testim = zeros(size(E,1),size(E,2));
%         testim(sub2ind(size(E),points(inliers,2),points(inliers,1))) = 1;
%         figure;imshow(testim);
%         end 
        
        p_dot_temp = dot(normals(inliers,:), ellipse_normals, 2);
        p_cnt = sum(p_dot_temp>0);%����֮�٣���һ�μ���ͳ�ƣ�����ΪC����ʱע�����Բʱ�ڵ㼫�Ե�ѡȡ����
        if(p_cnt > size(inliers,1)*0.5)
            %��������,Ҳ�����ں����    
            ellipse_polarity = -1;
            inliers = inliers(p_dot_temp>0 & p_dot_temp >= 0.923879532511287 );%cos(pi/8) = 0.923879532511287, �н�С��22.5��  
        else
            %������ͬ,Ҳ�����ڰ���� 
            ellipse_polarity = 1;
            inliers = inliers(p_dot_temp<0 & (-p_dot_temp) >= 0.923879532511287 );
        end
        
%         if( i == 11 && angleCoverage == 165)
%         testim = zeros(size(E,1),size(E,2));
%         testim(sub2ind(size(E),points(inliers,2),points(inliers,1))) = 1;
%         figure;imshow(testim);
%         end
        
        inliers2 = inliers;
        inliers3 = 0;
        %=================================================================================================================================================================
        %��ͨ�������inliersΪ������ڱ�Ե������±�
%         size(points)
%          size(inliers)
%         size(points(inliers, :))
%         size(takeInliers(points(inliers, :), circleCenter, tbins))
        %��ͨ��������õ���Ч���ڵ㣬�ڵ��ᴿ��Ҳ����inliers�н�һ��������Ч��inliers,��������٣���Сinlier_n2 x 1��ע��ע�⣺inliers�д������points�е����±�
        inliers = inliers(takeInliers(points(inliers, :), ellipseCenter, tbins));
        
%         if( i == 11 && angleCoverage == 165)
%         testim = zeros(size(E,1),size(E,2));
%         testim(sub2ind(size(E),points(inliers,2),points(inliers,1))) = 1;
%         figure;imshow(testim);
%         end

         [new_ellipse,new_info] = fitEllipse(points(inliers,1),points(inliers,2));
         
%           if( i == 11 && angleCoverage == 165)
%          drawEllipses(new_ellipse',E);
%           end

%         if angleLoop == 2   %mycode
%         dispimg = zeros(size(E,1),size(E,2),3);
%         dispimg(:,:,1) = E.*255;%��Ե��ȡ��������0-1ͼ��
%         dispimg(:,:,2) = E.*255;
%         dispimg(:,:,3) = E.*255;
%         for i = 1:length(inliers)
%         dispimg(points(inliers(i),2),points(inliers(i),1),:)=[0 0 255];
%         end
%         dispimg = drawCircle(dispimg,[newa newb],newr);
%         figure;
%         imshow(uint8(dispimg));
%         end

        if (new_info == 1)%���������С���˷���ϵĶ��ó��Ľ��
            %�¶Գ����ĺ��϶Գ����ĵľ���С��4*distance_tolerance, (a,b)�ľ���Ҳ��С��4*distance_tolerance,���phiС��0.314159265358979 = 0.1pi = 18��,��Ϊ����ϳ����Ĳ��ܺ�ԭ������Բ���Ĳ�ܶ�,
            if ( (((new_ellipse(1) - ellipseCenter(1))^2 + (new_ellipse(2) - ellipseCenter(2))^2 ) <= 16 * distance_tolerance_square) ...
                && (((new_ellipse(3) - ellipseAxes(1))^2 + (new_ellipse(4) - ellipseAxes(2))^2 ) <= 16 * distance_tolerance_square) ...
                && (abs(new_ellipse(5) - ellipsePhi) <= 0.314159265358979) )
                ellipse_normals = computePointAngle(new_ellipse,points);
                %������һ�����ڵ㣬��ͨ�Է������ڵ��ᴿ,��ε��µ��ڵ�����ں���������ȷ���
                %newinliers = find( (labels == 0) & (dRosin_square(new_ellipse,points) <= distance_tolerance_square) ...
                %    & ((dot(normals, ellipse_normals, 2)*(-ellipse_polarity)) >= 0.923879532511287) ); % (2*distance_tolerance)^2, cos(pi/8) = 0.923879532511287, �н�С��22.5��
                %���ټ��㣬ֻ������Բ��Ӿ����ڵı�Ե��(��Բ�еĳ���a>b),i_dx�洢���������points�е�����
                i_dx = find( points(:,1) >= (new_ellipse(1)-new_ellipse(3)-distance_tolerance) & points(:,1) <= (new_ellipse(1)+new_ellipse(3)+distance_tolerance) & points(:,2) >= (new_ellipse(2)-new_ellipse(3)-distance_tolerance) & points(:,2) <= (new_ellipse(2)+new_ellipse(3)+distance_tolerance));
                ellipse_normals = computePointAngle(new_ellipse,points(i_dx,:));%ellipse_normals������i_dx����һ��
                newinliers = i_dx(labels(i_dx) == 0 & (dRosin_square(new_ellipse,points(i_dx,:)) <= distance_tolerance_square & ((dot(normals(i_dx,:), ellipse_normals, 2)*(-ellipse_polarity)) >= 0.923879532511287) ) );
                newinliers = newinliers(takeInliers(points(newinliers, :), new_ellipse(1:2), tbins));
                if (length(newinliers) >= length(inliers))
                    %a = newa; b = newb; r = newr; cnd = newcnd;
                    inliers = newinliers;
                    inliers3 = newinliers;%my code��just test
                    %======================================================================
                    %�������
                    %[newa, newb, newr, newcnd] = fitCircle(points(inliers, :));
                    [new_new_ellipse,new_new_info] = fitEllipse(points(inliers,1),points(inliers,2));
                    if(new_new_info == 1)
                       new_ellipse = new_new_ellipse;
                    end
                    %=======================================================================
                end
            end
        else
            new_ellipse = list(i,:);  %candidates
        end
        
        %�ڵ���������Բ���ϵ�һ��������TminΪ������ֵ
%         length(inliers)
%         floor( pi * (1.5*sum(new_ellipse(3:4))-sqrt(new_ellipse(3)*new_ellipse(4))) * Tmin )
        if (length(inliers) >= floor( pi * (1.5*sum(new_ellipse(3:4))-sqrt(new_ellipse(3)*new_ellipse(4))) * Tmin ))
            convergence(i, :) = new_ellipse;
            %��֮ǰ��Բ�ĺͰ뾶��������һ�£��ظ��ˣ���˰����Բ��̭(�����ͷ�ĺ����ظ���Բ��һ���ᱻ��̭)
            if (any( (sqrt(sum((convergence(1 : i - 1, 1 : 2) - repmat(new_ellipse(1:2), i - 1, 1)) .^ 2, 2)) <= distance_tolerance) ...
                & (sqrt(sum((convergence(1 : i - 1, 3 : 4) - repmat(new_ellipse(3:4), i - 1, 1)) .^ 2, 2)) <= distance_tolerance) ...
                & (abs(convergence(1 : i - 1, 5) - repmat(new_ellipse(5), i - 1, 1)) <= 0.314159265358979) ))
                validCandidates(i) = false;
            end
            %����ڵ���Բ��������angleCoverage��������
            %completeOrNot =  isComplete(points(inliers, :), new_ellipse(1:2), tbins, angleCoverage);
            completeOrNot = calcuCompleteness(points(inliers,:),new_ellipse(1:2),tbins) >= angleCoverage;
            if (new_info == 1 && new_ellipse(3) < maxSemiMajor && new_ellipse(4) < maxSemiMinor && completeOrNot )
                %�����������Բ��������distance_tolerance��Ҳ����ָ������Բ�ǲ�ͬ��
                if (all( (sqrt(sum((ellipses(:, 1 : 2) - repmat(new_ellipse(1:2), size(ellipses, 1), 1)) .^ 2, 2)) > distance_tolerance) ...
                   | (sqrt(sum((ellipses(:, 3 : 4) - repmat(new_ellipse(3:4), size(ellipses, 1), 1)) .^ 2, 2)) > distance_tolerance) ...
                   | (abs(ellipses(:, 5) - repmat(new_ellipse(5), size(ellipses, 1), 1)) >= 0.314159265358979 ) )) %0.1 * pi = 0.314159265358979 = 18��
                    %size(inliers)
                    %line_normal = pca(points(inliers, :));%�õ�2x2��pca�任������˵ڶ��б������ڵ�ͳ�Ƴ����ݶ�
                    %line_normal = line_normal(:, 2)';%ȡ���ڶ��в��ұ�Ϊ1 x 2 ��������
                    %line_point = mean(points(inliers, :));%�ڵ�ȡƽ��
                    %��ֹ���ݵ���ڼ���
                    %if (sum(abs(dot(points(inliers, :) - repmat(line_point, length(inliers), 1), repmat(line_normal, length(inliers), 1), 2)) <= distance_tolerance & radtodeg(real(acos(abs(dot(normals(inliers, :), repmat(line_normal, length(inliers), 1), 2))))) <= normal_tolerance) / length(inliers) < 0.8)
                         labels(inliers) = size(ellipses, 1) + 1;%��ǣ���Щ�ڵ��Ѿ��ù��ˣ��������¼�⵽Բ��
                         %==================================================================
                         if(all(inliers3) == 1)
                         mylabels(inliers3) = size(ellipses,1) + 1; %��ʾ��ϳ�Բʱ���õ��ڵ�  SSS
                         end
                         %==================================================================
                        ellipses = [ellipses; new_ellipse];%����Բ���������ȥ
                        validCandidates(i) = false;%��i����ѡԲ������
                        %disp([angleCoverage,i]);
                        %drawEllipses(new_ellipse',E);
                    %end
                end
            end
        else
            validCandidates(i) = false;%�����������̭�ú�ѡԲ
        end
        
    end %for
end%fun
%% ================================================================================================================================
%����4
%Բ����С���˷����(�˴����Ը��ÿ���Բ��Ϸ���)
%���룺
%points: ��ͨ�Է�������ᴿ����ڵ㣬���СΪ fpn x 2,��ʽ(xi,yi)
%�����
%a   ����Ϻ��Բ�ĺ�����x
%b   ����Ϻ��Բ��������y
%c   ����Ϻ��Բ�İ뾶r
%cnd ��1��ʾ���ݴ��뷽�̺�������ģ�ֱ����ƽ��ֵ���ƣ�0��ʾ����������С���˷���ϵ�
function [a, b, r, cnd] = fitCircle(points)
%{
    A = [sum(points(:, 1)), sum(points(:, 2)), size(points, 1); sum(points(:, 1) .* points(:, 2)), sum(points(:, 2) .* points(:, 2)), sum(points(:, 2)); sum(points(:, 1) .* points(:, 1)), sum(points(:, 1) .* points(:, 2)), sum(points(:, 1))];
    %����С���˷�ʱ��A'A�����������ӽ�0������ζ�ŷ��������ԣ���ƽ��ֵ����
    if (abs(det(A)) < 1e-9)
        cnd = 1;
        a = mean(points(:, 1));
        b = mean(points(:, 2));
        r = min(max(points) - min(points));
        return;
    end
    cnd = 0;
    B = [-sum(points(:, 1) .* points(:, 1) + points(:, 2) .* points(:, 2)); -sum(points(:, 1) .* points(:, 1) .* points(:, 2) + points(:, 2) .* points(:, 2) .* points(:, 2)); -sum(points(:, 1) .* points(:, 1) .* points(:, 1) + points(:, 1) .* points(:, 2) .* points(:, 2))];
    t = A \ B;
    a = -0.5 * t(1);
    b = -0.5 * t(2);
    r = sqrt((t(1) .^ 2 + t(2) .^ 2) / 4 - t(3));
 %}
    A = [sum(points(:, 1) .* points(:, 1)),sum(points(:, 1) .* points(:, 2)),sum(points(:, 1)); sum(points(:, 1) .* points(:, 2)),sum(points(:, 2) .* points(:, 2)),sum(points(:, 2)); sum(points(:, 1)),sum(points(:, 2)),size(points, 1)]; 
    %����С���˷�ʱ��A'A�����������ӽ�0������ζ�ŷ��������ԣ���ƽ��ֵ����
    if (abs(det(A)) < 1e-9)
        cnd = 1;
        a = mean(points(:, 1));
        b = mean(points(:, 2));
        r = min(max(points) - min(points));
        return;
    end
    cnd = 0;
    B = [sum(-points(:, 1) .* points(:, 1) .* points(:, 1) - points(:, 1) .* points(:, 2) .* points(:, 2));sum(-points(:, 1) .* points(:, 1) .* points(:, 2) - points(:, 2) .* points(:, 2) .* points(:, 2)); sum(-points(:, 1) .* points(:, 1) - points(:, 2) .* points(:, 2))];
    t = A \ B;
    a = -0.5 * t(1);
    b = -0.5 * t(2);
    r = sqrt((t(1) .^ 2 + t(2) .^ 2) / 4 - t(3));
end
%% ================================================================================================================================
%����5
%����
%x     : ��ͨ�Է�������������2piRT���ᴿ����ڵ�(x,y)�������뵽�����ȷ�������.num x 2
%center: Բ��(x,y)  1 x 2
%tbins ����������
%angleCoverage: ��Ҫ�ﵽ��Բ������
%���
%result�� true or false����ʾ��Բ�����벻����
%longest_inliers:
function [result, longest_inliers] = isComplete(x, center, tbins, angleCoverage)
    [theta, ~] = cart2pol(x(:, 1) - center(1), x(:, 2) - center(2));%thetaΪ(-pi,pi)�ĽǶȣ�num x 1
    tmin = -pi; tmax = pi;
    tt = round((theta - tmin) / (tmax - tmin) * tbins + 0.5);%theta�ĵ�i��Ԫ�����ڵ�j��bin����tt��i�б��Ϊj����Сnum x 1
    tt(tt < 1) = 1; tt(tt > tbins) = tbins;
    h = histc(tt, 1 : tbins);
    longest_run = 0;
    start_idx = 1;
    end_idx = 1;
    while (start_idx <= tbins)
        if (h(start_idx) > 0)%�ҵ�bin��vote��һ������0��
            end_idx = start_idx;
            while (start_idx <= tbins && h(start_idx) > 0)%ֱ��bin��һ��С��0��
                start_idx = start_idx + 1;
            end
            inliers = [end_idx, start_idx - 1];%������Ϊ��ͨ����
            inliers = find(tt >= inliers(1) & tt <= inliers(2));%��tt���ҵ����ڴ�������ڵ������
            run = max(theta(inliers)) - min(theta(inliers));%�ǶȲ�
            if (longest_run < run)%�˾���Ϊ���ҵ���������������ͨ�Ŀ��
                longest_run = run;
                longest_inliers = inliers;
            end
        end
        start_idx = start_idx + 1;
    end
    if (h(1) > 0 && h(tbins) > 0)%�����һ��bin�����һ��bin������0���п��������ͨ������ͷβ�������������
        start_idx = 1;
        while (start_idx < tbins && h(start_idx) > 0)%�ҵ�bin��vote��һ������0��
            start_idx = start_idx + 1;
        end
        end_idx = tbins;%end_idxֱ�Ӵ���β����ʼ������
        while (end_idx > 1 && end_idx > start_idx && h(end_idx) > 0)
            end_idx = end_idx - 1;
        end
        inliers = [start_idx - 1, end_idx + 1];
        run = max(theta(tt <= inliers(1)) + 2 * pi) - min(theta(tt >= inliers(2)));
        inliers = find(tt <= inliers(1) | tt >= inliers(2));
        if (longest_run < run)
            longest_run = run;
            longest_inliers = inliers;
        end
    end
    %������ͨ�Ŀ�ȴ�����angleCoverage��������Ȼ�����ͨ���С�ڣ����������㹻��
    longest_run_deg = rad2deg(longest_run);
    h_greatthanzero_num = sum(h>0);
    result =  longest_run_deg >= angleCoverage || h_greatthanzero_num * (360 / tbins) >= min([360, 1.2*angleCoverage]);  %1.2 * angleCoverage
end
function [completeness] = calcuCompleteness(x, center, tbins)
    [theta, ~] = cart2pol(x(:, 1) - center(1), x(:, 2) - center(2));%thetaΪ(-pi,pi)�ĽǶȣ�num x 1
    tmin = -pi; tmax = pi;
    tt = round((theta - tmin) / (tmax - tmin) * tbins + 0.5);%theta�ĵ�i��Ԫ�����ڵ�j��bin����tt��i�б��Ϊj����Сnum x 1
    tt(tt < 1) = 1; tt(tt > tbins) = tbins;
    h = histc(tt, 1 : tbins);
    h_greatthanzero_num = sum(h>0);
    completeness = h_greatthanzero_num*(360 / tbins);
end
%% ================================================================================================================================
%����6
%��ͨ�Է�������Բ���ϵ��ڵ�����ᴿ
%����
%x����Բ���ϵ��ڵ�(x,y),��Ϊinlier_n x 2 
%center��һ����Բ������(x,y) 1x2
%tbins: ���� = min( 180 , pi*(1.5*(a+b)-sqrt(a*b)) ) 
%���
%idx��Ϊ��xһ�����ģ�inlier_n x 1��logical������������Ч������һ����ͨ���ȵ��ڵ㣬��Ӧλ����Ч��Ϊ1������Ϊ0
function idx = takeInliers(x, center, tbins)
   [theta, ~] = cart2pol(x(:, 1) - center(1), x(:, 2) - center(2));%�õ�[-pi,pi]�ķ�λ�ǣ��ȼ��� theta = atan2(x(:, 2) - center(2) , x(:, 1) - center(1)); 
    tmin = -pi; tmax = pi;
    tt = round((theta - tmin) / (tmax - tmin) * tbins + 0.5);%���ڵ������[1 tbins]
    tt(tt < 1) = 1; tt(tt > tbins) = tbins;
    h = histc(tt, 1 : tbins);%hΪֱ��ͼ[1 tbins]��ͳ�ƽ��
    mark = zeros(tbins, 1);
    compSize = zeros(tbins, 1);
    nComps = 0;
    queue = zeros(tbins, 1);
    du = [-1, 1];
    for i = 1 : tbins
        if (h(i) > 0 && mark(i) == 0)%������ڵ�i�������ڵ�ֵ����0����mark(i)Ϊ0
            nComps = nComps + 1;
            mark(i) = nComps;%��ǵ�nComps����ͨ����
            front = 1; rear = 1;
            queue(front) = i;%���÷���������У����Դ˿�ʼ����
            while (front <= rear)
                u = queue(front);
                front = front + 1;
                for j = 1 : 2
                    v = u + du(j);
                    if (v == 0)
                        v = tbins;
                    end
                    if (v > tbins)
                        v = 1;
                    end
                    if (mark(v) == 0 && h(v) > 0)
                        rear = rear + 1;
                        queue(rear) = v;
                        mark(v) = nComps;%��ǵ�nComps����ͨ����
                    end
                end
            end
            compSize(nComps) = sum(ismember(tt, find(mark == nComps)));%�õ�������ͨ��ΪnComps���ڵ�����
        end
    end
    compSize(nComps + 1 : end) = [];
    maxCompSize = max(compSize);
    validComps = find(compSize >= maxCompSize * 0.1 & compSize > 10);%���ڵ��������ͨ���ȵ�0.1������ͨ��������Ч��
    validBins = find(ismember(mark, validComps));%��Ч�ķ���
    idx = ismember(tt, validBins);%��Ч���ڵ�
end
%% compute the points' normals belong to an ellipse, the normals have been already normalized. 
%param: [x0 y0 a b phi].
%points: [xi yi], n x 2
function [ellipse_normals] = computePointAngle(ellipse, points)
%convert [x0 y0 a b phi] to Ax^2+Bxy+Cy^2+Dx+Ey+F = 0
a_square = ellipse(3)^2;
b_square = ellipse(4)^2;
sin_phi = sin(ellipse(5));
cos_phi = cos(ellipse(5)); 
sin_square = sin_phi^2;
cos_square = cos_phi^2;
A = b_square*cos_square+a_square*sin_square;
B = (b_square-a_square)*sin_phi*cos_phi*2;
C = b_square*sin_square+a_square*cos_square;
D = -2*A*ellipse(1)-B*ellipse(2);
E = -2*C*ellipse(2)-B*ellipse(1);
% F = A*ellipse(1)^2+C*ellipse(2)^2+B*ellipse(1)*ellipse(2)-(ellipse(3)*ellipse(4)).^2;
% A = A/F;
% B = B/F;
% C = C/F;
% D = D/F;
% E = E/F;
% F = 1;
%calculate points' normals to ellipse
angles = atan2(C*points(:,2)+B/2*points(:,1)+E/2, A*points(:,1)+B/2*points(:,2)+D/2);
ellipse_normals = [cos(angles),sin(angles)];
end

%% paramΪ[x0 y0 a b Phi],1 x 5 ���� 5 x 1
%pointsΪ������rosin distance�ĵ㣬ÿһ��Ϊ(xi,yi),size�� n x 2
%dminΪ������ƾ����ƽ��.
%����ע�⣬��a = bʱ��Ҳ������Բ�˻���Բʱ��dmin���������NAN�������ô˺���
function [dmin]= dRosin_square(param,points)
ae2 = param(3).*param(3);
be2 = param(4).*param(4);
x = points(:,1) - param(1);
y = points(:,2) - param(2);
xp = x*cos(-param(5))-y*sin(-param(5));
yp = x*sin(-param(5))+y*cos(-param(5));
fe2 = ae2-be2;
X = xp.*xp;
Y = yp.*yp;
delta = (X+Y+fe2).^2-4*fe2*X;
A = (X+Y+fe2-sqrt(delta))/2;
ah = sqrt(A);
bh2 = fe2-A;
term = A*be2+ae2*bh2;
xi = ah.*sqrt(ae2*(be2+bh2)./term);
yi = param(4)*sqrt(bh2.*(ae2-A)./term);
d = zeros(size(points,1),4);%n x 4
d(:,1) = (xp-xi).^2+(yp-yi).^2;
d(:,2) = (xp-xi).^2+(yp+yi).^2;
d(:,3) = (xp+xi).^2+(yp-yi).^2;
d(:,4) = (xp+xi).^2+(yp+yi).^2;
dmin = min(d,[],2); %���ؾ����ƽ��
%[dmin, ii] = min(d,[],2); %���ؾ����ƽ��
% for jj = 1:length(dmin)
%     if(ii(jj) == 1)
%         xi(jj) = xi(jj);
%         yi(jj) = yi(jj);
%     elseif (ii(jj) == 2)
%         xi(jj) = xi(jj);
%         yi(jj) = -yi(jj);
%     elseif (ii(jj) == 3)
%         xi(jj) = -xi(jj);
%         yi(jj) = yi(jj);
%     elseif(ii(jj) == 4)
%          xi(jj) = -xi(jj);
%         yi(jj) = -yi(jj);
%     end
% end
% 
% xi =  xi*cos(param(5))-yi*sin(param(5));
% yi =  xi*sin(param(5))+yi*cos(param(5));
% 
% testim = zeros(300,300);
% testim(sub2ind([300 300],uint16(yi+param(2)),uint16(xi+param(1)))) = 1;
% figure;imshow(uint8(testim).*255);
end
