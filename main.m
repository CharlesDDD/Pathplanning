clear;
close all;
tic%%时间记录
G = [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0
    0 1 1 0 1 0 0 0 0 0 1 0 0 0 0 1 1 1 0 0
    0 1 1 0 1 0 1 1 1 0 0 1 0 0 0 1 1 1 0 0
    0 1 0 0 0 1 1 1 1 1 0 0 1 0 0 1 1 1 0 0
    0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 1 1 1 0 0
    0 0 0 1 0 0 1 1 1 0 0 1 0 0 0 1 1 1 0 0
    0 1 1 1 0 0 1 1 1 0 0 0 0 1 0 1 1 1 0 0
    0 1 1 1 1 0 1 1 1 0 1 1 1 1 0 1 1 1 0 0
    0 1 1 1 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0
    0 1 1 1 0 0 0 0 0 0 1 1 1 1 0 0 1 0 0 0
    0 0 0 0 1 0 0 1 1 0 1 1 1 1 1 0 0 0 0 0
    0 0 0 1 1 1 0 1 1 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 1 0 0 0 0 0 1 0 1 1 1 0 1 1 1 1 0
    0 0 0 1 0 0 0 0 0 1 0 1 1 1 0 1 1 1 1 0
    0 0 0 0 0 0 0 1 0 0 0 1 1 1 0 1 1 1 1 0
    0 0 1 1 0 0 1 1 1 0 0 0 0 1 1 0 0 0 0 0
    0 0 1 1 0 0 1 1 1 0 1 1 0 0 0 0 1 1 1 0
    0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 0 0 1 1 0
    0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 1 0;];


MM=size(G,1);% G 地形图为01矩阵，如果为1表⽰障碍物
Tau=ones(MM*MM,MM*MM);% Tau 初始信息素矩阵
% Tau=8.*Tau;
K=100;%迭代次数
M=50;%蚂蚁个数
S=1;%最短路径的起始点 此时的起点在左下角坐标
E=400;%最短路径的⽬的点 从左到右从上到下 此时的终点在右上角此时的坐标
% q0 = 0.5;
%%[row1, col1][row2, col2]这里求出来起止点的索引ind对应的行列 这里是按照矩阵规定的上下左右来的，用来让起止点为自由点
sz = [MM MM];
[row1, col1] = ind2sub(sz, S);
[row2, col2] = ind2sub(sz, E);

%%随机生成障碍物
% for i = 1:MM
%     for j = 1:MM
%         G(i,j) = 1-(rand(1)>=.5)*1;
%         if (i==row1&&j==col1) || (i==row2&&j==col2)%%起止点一定要是0可行区域否则会偶尔有bug
%             G(i,j) = 0;
%         end
%         if G(i,j) == 1
%             G(i,j) = 1-(rand(1)>=.5)*1;
%         end
%     end
% end

Alpha=0.9;% Alpha 表征信息素重要程度的参数
Beta=9;% Beta 表征启发式因⼦重要程度的参数
Rho=0.2;% Rho 信息素蒸发系数
Q=1;% Q 信息素增加强度系数
minkl=inf;
mink=0;
minl=0;
D=G2D(G);%%每个自由格栅到它48个邻域中的自由格栅
N=size(D,1);%N表⽰问题的规模（象素个数）
a=1;%⼩⽅格象素的边长

%%Ex Ey Sx Sy求出来是为了画起点终点的圆圈
Ex=a*(mod(E,MM)-0.5);%终点横坐标
if Ex==-0.5
    Ex=MM-0.5;
end
Ey=a*(MM+0.5-ceil(E/MM));%终点纵坐标
Sx=a*(mod(S,MM)-0.5);%起点横坐标
if Sx==-0.5
    Sx=MM-0.5;
end
Sy=a*(MM+0.5-ceil(S/MM));%起点纵坐标
D_SE = ((Sx-Ex)^2+(Sy-Ey)^2)^0.5;%%起止点距离
Eta=zeros(N);%启发式信息，取为当前格栅（从左到右从上到下）⾄⽬标点的直线距离的倒数--每一个格栅到终点格栅的距离的倒数
%以下启发式信息矩阵  
for i=1:N
    ix=a*(mod(i,MM)-0.5);
    if ix==-0.5
        ix=MM-0.5;
    end
    iy=a*(MM+0.5-ceil(i/MM));
    dist_SE = ((Sx-Ex)^2+(Sy-Ey)^2)^0.5;
    dist_Si = ((ix-Sx)^2+(iy-Sy)^2)^0.5;
    dist_iE = ((ix-Ex)^2+(iy-Ey)^2)^0.5;
    Tau(S,i) = dist_SE/(dist_Si + dist_iE);
    Tau(i,S) = dist_SE/(dist_Si + dist_iE);
    Tau(E,i) = dist_SE/(dist_Si + dist_iE);
    Tau(i,E) = dist_SE/(dist_Si + dist_iE);
    if i~=E
        Eta(i)=1/((ix-Ex)^2+(iy-Ey)^2)^0.5;
    else
        Eta(i)=100;
    end
end
%%不均匀分配信息素
for i = 1:N
    iix=a*(mod(i,MM)-0.5);
    if iix==-0.5
        iix=MM-0.5;
    end
    iiy=a*(MM+0.5-ceil(i/MM));
    for j = 1:N
        jjx=a*(mod(j,MM)-0.5);
        if jjx==-0.5
            jjx=MM-0.5;
        end
        jjy=a*(MM+0.5-ceil(j/MM));
        Tau(i,j) = 1/(((iix-jjx)^2+(iiy-jjy)^2)^0.5+((jjx-Ex)^2+(jjy-Ey)^2)^0.5)+Tau(i,j);
        if i == j
            Tau(i,j) = 0;
        end
    end
end
% 障碍物处的信息素设置为-1
for i= 1:MM
    for j = 1:MM
        if G(j,i) == 1
            ind = 20*(j-1) + i;
            Tau(ind,S) = -1;
            Tau(S,ind) = -1;
            Tau(ind,E) = -1;
            Tau(E,ind) = -1;
        end
    end
end

ROUTES=cell(K,M);%⽤细胞结构存储每⼀代的每⼀只蚂蚁的爬⾏路线
PL=zeros(K,M);%⽤矩阵存储每⼀代的每⼀只蚂蚁的爬⾏路线长度
%启动K轮蚂蚁觅⾷活动，每轮派出M只蚂蚁
for k=1:K
    TauStart=Tau;
    for m=1:M
        %状态初始化
        W=S;%当前节点初始化为起始点
        Path=S;%爬⾏路线初始化
        PLkm=0;%爬⾏路线长度初始化
        TABUkm=ones(N);%禁忌表初始化 没有走过的为1 走过的设置为0  TABUkm只用到了第一列
        TABUkm(S)=0;%已经在初始点了，因此要排除
        DD=D;%邻接矩阵初始化
        %下⼀步可以前往的节点
        DW=DD(W,:);
        DW1=find(DW);%%找到非0的索引值
        for j=1:length(DW1)
            if TABUkm(DW1(j))==0
                DW(DW1(j))=0;
            end
        end
        LJD=find(DW);
        Len_LJD=length(LJD);%可选节点的个数
        %蚂蚁未遇到⾷物或者陷⼊死胡同或者觅⾷停⽌
        while W~=E&&Len_LJD>=1 %%W~=E没到终点 Len_LJD>=1至少有一个可选节点
            %转轮赌法选择下⼀步怎么⾛
            PP=zeros(Len_LJD);
            for i=1:Len_LJD
                %%Tau(W,LJD(i))表示当前节点到下一个可走节点的信息素浓度Eta(LJD(i)是下一个节点到终点距离的倒数
                %%传统Eta
                PP(i)=(Tau(W,LJD(i))^Alpha)*((Eta(LJD(i)))^Beta);

                %%改进...
                %%当前节点W
%                 if K>1
%                     Wx=a*(mod(W,MM)-0.5);
%                     if Wx==-0.5
%                         Wx=MM-0.5;
%                     end
%                     Wy=a*(MM+0.5-ceil(W/MM));
%                     %%当前节点下一个节点LJD(i)
%                     temp = LJD(i);
%                     temp_x=a*(mod(temp,MM)-0.5);
%                     if temp_x==-0.5
%                         temp_x=MM-0.5;
%                     end
%                     temp_y=a*(MM+0.5-ceil(temp/MM));
%                     D_SW = ((Sx-Wx)^2+(Sy-Wy)^2)^0.5;
%                     D_Wtwmp = ((Wx-temp_x)^2+(Wy-temp_y)^2)^0.5;
%                     D_tempE = ((temp_x-Ex)^2+(temp_y-Ey)^2)^0.5;
%                     PP(i)=(Tau(W,LJD(i))^Alpha)*((1/(D_SW+D_Wtwmp+D_tempE))^Beta);
%                 else
%                     PP(i)=(Tau(W,LJD(i))^Alpha)*((Eta(LJD(i)))^Beta);
%                 end
            end
            sumpp=sum(PP);
            PP=PP/sumpp;%建⽴概率分布

            %%改进伪随机概率2
            if k>= 2 && k <=60
                 q = rand;
                 q0 = min(nonzeros(PL(K-1,:)))/max(PL(K-1,:));
                if q <= q0
                    Select=find(PP==max(PP));
                else
                    Pcum = cumsum(PP);
                    Select=find(Pcum>=rand);
                end
            elseif k<2
                Pcum = cumsum(PP);
                Select=find(Pcum>=rand);
            else
                q = rand;
                q0 = D_SE/minkl;
                if q <= q0
                    Select=find(PP==max(PP));
                else
                    Pcum = cumsum(PP);
                    Select=find(Pcum>=rand);
                end
            end
           
            %%传统概率
%             Pcum = cumsum(PP);
%             Select=find(Pcum>=rand);

            %%传统伪随机概率
%             q = rand;
%             q0 = min(nonzeros(PL(K-1,:)))/max(PL(K-1,:));
%             if q <= q0
%                 Select=find(PP==max(PP));
%             else
%                 Pcum = cumsum(PP);
%                 Select=find(Pcum>=rand);
%             end

            to_visit=LJD(Select(1));
            %状态更新和记录
            Path=[Path,to_visit];%路径增加
            PLkm=PLkm+DD(W,to_visit);%路径长度增加
            W=to_visit;%蚂蚁移到下⼀个节点
            %%这里移动过之后不往回退，走过的代价变为0，
            for kk=1:N
                if TABUkm(kk)==0
                    DD(W,kk)=0;
                    DD(kk,W)=0;%%设置为0，避免回退
                end
            end
            TABUkm(W)=0;%已访问过的节点从禁忌表中删除
            DW=DD(W,:);
            DW1=find(DW);
            for j=1:length(DW1)
                if TABUkm(DW1(j))==0 %%如果走过了
                    DW(j)=0;%%就从当前点的可走邻域中去除，只要没走过的邻域
                end
            end
            LJD=find(DW);
            Len_LJD=length(LJD);%可选节点的个数
        end
        %记下每⼀代每⼀只蚂蚁的觅⾷路线和路线长度
        ROUTES{k,m}=Path;
        if Path(end)==E%判断本只蚂蚁寻找路径的最后一个节点是否为终点
            PL(k,m)=PLkm;
            if PLkm<minkl
                mink=k;minl=m;minkl=PLkm;%%最短路径是在那一代的哪一只蚂蚁
            end
        else
%             LastNode = W;
%             Path_size = size(Path, 2);
%             Rand_path = randi([1, Path_size-1], 1, 1);
%             Start_node = Path(1, Rand_path);
%             W = Start_node;
%             Path = Path(1:Rand_path);
%             PLkm = sum(DD(Path(1:end-1), Path(2:end)));
%             
%             if Path(end) == E
%                 PL(k,m) = PLkm;
%                 if PLkm < minkl
%                     mink = k; minl = m; minkl = PLkm;
%                 end
%             else
%                 W = LastNode;
%                 PL(k,m)=0;
%             end
            %%二次路径规划
            random_index = randi(length(Path) - 1); 
            W = Path(random_index); 
            Path = Path(1:random_index);
            DW = DD(W,:);
            DW1 = find(DW);
            for j = 1:length(DW1)
                if TABUkm(DW1(j)) == 0
                    DW(DW1(j)) = 0;
                end
            end
            LJD = find(DW);
            Len_LJD = length(LJD);
            while W~=E&&Len_LJD>=1 %%W~=E没到终点 Len_LJD>=1至少有一个可选节点
                %转轮赌法选择下⼀步怎么⾛
                PP=zeros(Len_LJD);
                for i=1:Len_LJD
                    %%Tau(W,LJD(i))表示当前节点到下一个可走节点的信息素浓度Eta(LJD(i)是下一个节点到终点距离的倒数
                    %%传统Eta
                    PP(i)=(Tau(W,LJD(i))^Alpha)*((Eta(LJD(i)))^Beta);
                end
                sumpp=sum(PP);
                PP=PP/sumpp;

                %%传统概率
                Pcum = cumsum(PP);
                Select=find(Pcum>=rand);

                to_visit=LJD(Select(1));
                %状态更新和记录
                Path=[Path,to_visit];
                PLkm=PLkm+DD(W,to_visit);
                W=to_visit;
                %%这里移动过之后不往回退，走过的代价变为0，
                for kk=1:N
                    if TABUkm(kk)==0
                        DD(W,kk)=0;
                        DD(kk,W)=0;
                    end
                end
                TABUkm(W)=0;%已访问过的节点从禁忌表中删除
                DW=DD(W,:);
                DW1=find(DW);
                for j=1:length(DW1)
                    if TABUkm(DW1(j))==0
                        DW(j)=0;
                    end
                end
                LJD=find(DW);
                Len_LJD=length(LJD);%可选节点的个数
            end
            %记下每⼀代每⼀只蚂蚁的觅⾷路线和路线长度
            ROUTES{k,m}=Path;
            if Path(end)==E%判断本只蚂蚁寻找路径的最后一个节点是否为终点
                PL(k,m)=PLkm;
                if PLkm<minkl
                    mink=k;minl=m;minkl=PLkm;%%最短路径是在那一代的哪一只蚂蚁
                end
            else
                 PL(k,m)=0;
            end
%             PL(k,m)=0;
        end
    end
    %%记录每代蚂蚁的死亡数
    c = 0;
    for i = 1:K
        for j = 1:M
            if PL(i,j) == 0
                c = c+1;
            end
        end
        death(i) = c;
        c= 0;
    end
    min_death = sum(death);
    %更新信息素
    Delta_Tau=zeros(N,N);%更新量初始化
    for m=1:M
        if PL(k,m)%%只更新到终点的蚂蚁的路径
            ROUT=ROUTES{k,m};
            TS=length(ROUT)-1;%跳数即：每两点之间是一个跳数
            PL_km=PL(k,m);%%当前蚂蚁走到终点的路径作为Q下面的分母
            for s=1:TS
                x=ROUT(s);
                y=ROUT(s+1);
%                 Delta_Tau(x,y)=Delta_Tau(x,y)+Q/PL_km + (Q*min(PL))/max(PL);
%                 Delta_Tau(y,x)=Delta_Tau(y,x)+Q/PL_km + (Q*min(PL))/max(PL);
                Delta_Tau(x,y)=Delta_Tau(x,y)+Q/PL_km;
                Delta_Tau(y,x)=Delta_Tau(y,x)+Q/PL_km;
            end
        end
    end
    Tau=(1-Rho).*Tau+Delta_Tau;%信息素挥发⼀部分，新增加⼀部分
    %%最大最小蚂蚁系统
    temp_1 = Tau; 
    Tau_max = M/(Rho*minkl);
    Tau_min = -(Tau_max/(0.1*M));
    for Row = 1:N
        for col = 1:N
            if Tau(Row,col) > Tau_max
                Tau(Row,col) = Tau_max;
            elseif Tau(Row,col) < Tau_min
                Tau(Row,col) = Tau_min;
            else
                Tau(Row,col) = temp_1(Row,col);
            end
        end
    end    
end


time = toc;%%时间记录
%绘图
plotif=1;%是否绘图的控制参数
if plotif==1%绘收敛曲线
%     minPL=zeros(K,1);
%     for i=1:K
%         PLK=PL(i,:);%%第一代50只蚂蚁各自走过的路径长度（只有到终点的）
%         Nonzero=find(PLK);%%找到到终点的蚂蚁
%         PLKPLK=PLK(Nonzero);%%把到终点的蚂蚁走过的路径长度存起来
%         minPL(i)=min(PLKPLK);%%选出本次迭代的到终点的蚂蚁的最短路径
%     end
%     figure(1)
%     plot(minPL);
%     hold on
%     grid on
% %     title('收敛曲线变化趋势');
%     xlabel('Number of iterations');
%     ylabel('The optimal path length');
    %路径绘制
    figure(2)
    axis([0,MM,0,MM])
    for i=1:MM
        for j=1:MM
            if G(i,j)==1
                x1=j-1;y1=MM-i;
                x2=j;y2=MM-i;
                x3=j;y3=MM-i+1;
                x4=j-1;y4=MM-i+1;
                fill([x1,x2,x3,x4],[y1,y2,y3,y4],[0.2,0.2,0.2]);
                hold on
            else
                x1=j-1;y1=MM-i;
                x2=j;y2=MM-i;
                x3=j;y3=MM-i+1;
                x4=j-1;y4=MM-i+1;
                fill([x1,x2,x3,x4],[y1,y2,y3,y4],[1,1,1]);
                hold on
            end
        end
    end
    hold on
    ROUT=ROUTES{mink,minl}; 
    LENROUT=length(ROUT);
    Rx=ROUT;
    Ry=ROUT;
    for ii=1:LENROUT
        Rx(ii)=a*(mod(ROUT(ii),MM)-0.5);
        if Rx(ii)==-0.5
            Rx(ii)=MM-0.5;
        end
        Ry(ii)=a*(MM+0.5-ceil(ROUT(ii)/MM));
    end
    plot(Rx,Ry, '-r', 'LineWidth', 2); %最短路径
    scatter(Sx, Sy, 50, 'filled', MarkerEdgeColor='b', MarkerFaceColor='red'); % 起点
    scatter(Ex, Ey, 50, 'filled', MarkerEdgeColor='r', MarkerFaceColor='blue'); % 终点

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%实验：记录第一代前10只蚂蚁走过的路径%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure(3)
%     axis([0,MM,0,MM])
%     for i=1:MM
%         for j=1:MM
%             if G(i,j)==1
%                 x1=j-1;y1=MM-i;
%                 x2=j;y2=MM-i;
%                 x3=j;y3=MM-i+1;
%                 x4=j-1;y4=MM-i+1;
%                 fill([x1,x2,x3,x4],[y1,y2,y3,y4],[0.2,0.2,0.2]);
%                 hold on
%             else
%                 x1=j-1;y1=MM-i;
%                 x2=j;y2=MM-i;
%                 x3=j;y3=MM-i+1;
%                 x4=j-1;y4=MM-i+1;
%                 fill([x1,x2,x3,x4],[y1,y2,y3,y4],[1,1,1]);
%                 hold on
%             end
%         end
%     end
%     n = 1;
%     for jj = 1:10
%         ROUT = ROUTES{1,jj};
%         len_rout = length(ROUT);
%         Rx=ROUT;
%         Ry=ROUT;
%         for ll=1:len_rout
%             Rx(ll)=a*(mod(ROUT(ll),MM)-0.5);
%             if Rx(ll)==-0.5
%                 Rx(ll)=MM-0.5;
%             end
%                 Ry(ll)=a*(MM+0.5-ceil(ROUT(ll)/MM));
%         end
%         switch n
%             case 1
%                 plot(Rx,Ry, '-g', 'LineWidth', 0.5);
%             case 2
%                 plot(Rx,Ry, '-b', 'LineWidth', 0.5);
%             case 3
%                 plot(Rx,Ry, '-r', 'LineWidth', 0.5);
%             case 4
%                 plot(Rx,Ry, '-c', 'LineWidth', 0.5);
%             case 5
%                 plot(Rx,Ry, '-m', 'LineWidth', 0.5);
%             case 6
%                 plot(Rx,Ry, '-y', 'LineWidth', 0.5);
%             case 7
%                 plot(Rx,Ry, 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 0.5);
%             case 8
%                 plot(Rx,Ry, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 0.5);
%             case 9
%                 plot(Rx,Ry, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 0.5);
%             case 10
%                 plot(Rx,Ry, 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 0.5);
%         end
%         n = n+1;
%           scatter(Sx, Sy, 50, 'filled', MarkerEdgeColor='b', MarkerFaceColor='red'); % 起点
%           scatter(Ex, Ey, 50, 'filled', MarkerEdgeColor='r', MarkerFaceColor='blue'); % 终点
%     end
end
% figure(4)
% grid on
% 
% plot(death,'r-');
% 
% %每代最差路径最优路径比值
figure(5)
grid on
for k = 1:60
    rat(k) = min(nonzeros(PL(k,:)))/max(PL(k,:));
end
for k = 61:100
    rat(k) = D_SE/min(nonzeros(PL(k,:)));
end
plot(rat,'bo-');


