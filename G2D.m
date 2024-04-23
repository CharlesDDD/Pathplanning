%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%计算各点的可行邻域%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D=G2D(G)
l=size(G,1);
D=zeros(l*l,l*l);%%D中的没20列代表G中的每一行
%%总共是400个点，计算每个点到每个点的代价
for i=1:l
    for j=1:l
        if G(i,j)==0%%前提是自由格栅
            for m=1:l
                for n=1:l
                    if G(m,n)==0%%计算的格栅也要是自由的
                        im=abs(i-m);
                        jn=abs(j-n);
                        %%判断是否是8邻域其中之一
                        if im+jn==1||(im==1&&jn==1)%%im+jn==1是8邻域上下左右 im==1&&jn==1是8邻域的4个对角
                            D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;%%
                        end

                        %%%%%%%%%%判断是否是16方向24邻域其中之一%%%%%%%%%%%%
                        %%1.四个直角方向
                        if im==2&&jn==2
                            if (m<i&&n<j) && G(m+1, n+1) == 0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;%%左上
                            elseif (m<i&&n>j) && G(m+1, n-1) == 0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;%%右上
                            elseif (m>i&&n<j) && G(m-1, n+1) == 0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;%%左下
                            elseif (m>i&&n>j) && G(m-1, n-1) == 0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;%%右下
                            end
                        end

                        %%直线四方向
                        if im+jn == 2
                            if (m<i&&n==j)&&G(m+1,n)==0%%上直
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            elseif (m>i&&n==j)&&G(m-1,n)==0%%下直
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            elseif (m==i&&n<j)&&G(m,n+1)==0%%左直
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            elseif (m==i&&n>j)&&G(m,n-1)==0%%右直
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            end
                        end

                        %%左右8对角的各两个2 3 6 7
                        if im==1&&jn==2
                            if i-m==1&&j-n==2&&G(m,n+1)==0&&G(i,j-1)==0%%左直上
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            elseif m-i==1&&j-n==2&&G(m,n+1)==0&&G(i,j-1)==0%%左直下
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            elseif m-i==1&&n-j==2&&G(m,n-1)==0&&G(i,j+1)==0%%右直下
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            elseif i-m==1&&n-j==2&&G(m,n-1)==0&&G(i,j+1)==0%%右直上
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            end
                        end

                        %%上下8对角的各两个1 4 5 8
                        if im==2&&jn==1
                            if i-m==2&&j-n==1&&G(m+1,n)==0&&G(i-1,j)==0%%上直左
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            elseif m-i==2&&j-n==1&&G(m-1,n)==0&&G(i+1,j)==0%%下直左
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            elseif m-i==2&&n-j==1&&G(m-1,n)==0&&G(i+1,j)==0%%下直右
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            elseif i-m==2&&n-j==1&&G(m+1,n)==0&&G(i-1,j)==0%%上直右
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            end
                        end

                        %%%%%%%%%%判断是否是32方向48邻域其中之一%%%%%%%%%%%%
                        %1.四个直角
                        if im==3&&jn==3
                            if (m<i&&n<j) && G(m+1, n+1) == 0 && G(i-1, j-1)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;%%左上
                            elseif (m<i&&n>j) && G(m+1, n-1) == 0 &&G(i-1, j+1)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;%%右上
                            elseif (m>i&&n<j) && G(m-1, n+1) == 0&&G(i+1, j-1)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;%%左下
                            elseif (m>i&&n>j) && G(m-1, n-1) == 0&&G(i+1, j+1)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;%%右下
                            end
                        end

                        %2.四个直线方向
                        if im+jn == 3
                            if (m<i&&n==j)&&G(m+1,n)==0 && G(i-1,j)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;%%上
                            elseif (m>i&&n==j)&&G(m-1,n)==0 && G(i+1,j)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;%%下
                            elseif (m==i&&n<j)&&G(m,n+1)==0 && G(i,j-1)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;%%左
                            elseif (m==i&&n>j)&&G(m,n-1)==0 && G(i,j+1)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;%%右
                            end
                        end

                        %3.斜角1
                        if im==2&&jn==3
                            if i-m==2&&j-n==3&&G(m,n+1)==0&&G(m+1,n+1)==0&&G(i,j-1)==0&&G(i-1,j-1)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            elseif m-i==2&&j-n==3&&G(m,n+1)==0&&G(m-1,n+1)==0&&G(i,j-1)==0&&G(i+1,j-1)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            elseif m-i==2&&n-j==3&&G(m,n-1)==0&&G(m-1,n-1)==0&&G(i,j+1)==0&&G(i+1,j+1)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            elseif i-m==2&&n-j==3&&G(m,n-1)==0&&G(m+1,n-1)==0&&G(i,j+1)==0&&G(i-1,j+1)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            end
                        end

                        %%4.斜角2
                        if im==3&&jn==2
                            if i-m==3&&j-n==2&&G(m+1,n)==0&&G(m+1,n+1)==0&&G(i-1,j)==0&&G(i-1,j-1)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            elseif m-i==3&&j-n==2&&G(m-1,n)==0&&G(m-1,n+1)==0&&G(i+1,j)==0&&G(i+1,j-1)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            elseif m-i==3&&n-j==2&&G(m-1,n)==0&&G(m-1,n-1)==0&&G(i+1,j)==0&&G(i+1,j+1)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            elseif i-m==3&&n-j==2&&G(m+1,n)==0&&G(m+1,n-1)==0&&G(i-1,j)==0&&G(i-1,j+1)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            end
                        end

                        %%5.斜角3
                        if im==1&&jn==3
                            if i-m==1&&j-n==3&&G(m,n+1)==0&&G(i,j-1)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            elseif m-i==1&&j-n==3&&G(m,n+1)==0&&G(i,j-1)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            elseif m-i==1&&n-j==3&&G(m,n-1)==0&&G(i,j+1)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            elseif i-m==1&&n-j==3&&G(m,n-1)==0&&G(i,j+1)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            end
                        end

                        %%6.斜角4
                        if im==3&&jn==1
                            if i-m==3&&j-n==1&&G(m+1,n)==0&&G(i-1,j)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            elseif m-i==3&&j-n==1&&G(m-1,n)==0&&G(i+1,j)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            elseif m-i==3&&n-j==1&&G(m-1,n)==0&&G(i+1,j)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            elseif i-m==3&&n-j==1&&G(m+1,n)==0&&G(i-1,j)==0
                                D((i-1)*l+j,(m-1)*l+n)=(im*im+jn*jn)^0.5;
                            end
                        end
                    end
                end
            end
        end
    end
end
end