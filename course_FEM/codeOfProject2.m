% A code calculate diplacement and stress using 2D 4-Node Isoparametric Formulation.
% It is for a FEM course in ZJU.

clc
clear;
node=[1 0 0 0;
      2 0.02 0 0;
      3 0.04 0 0;
      4 0 0.005 0;
      5 0.02 0.005 0;
      6 0.04 0.005 0;
      7 0 0.01 0;
      8 0.02 0.01 0;
      9 0.04 0.01 0];   %节点信息，第一列为节点编号，2~4列分别为x,y,z方向坐标
ele=[1 1 2 5 4;
     2 2 3 6 5;
     3 4 5 8 7;
     4 5 6 9 8];        %单元信息，第一列为单元编号，后面各列为单元上的节点号码
%---------------物理参数------------------
E=200e9;             %弹性模量
t=0.001;               %单元厚度
miu=0.3;           %泊松比
%-----------------------------------------
n_ele=length(ele(:,1));   %单元数

%组装总体刚度矩阵
dof=length(node(:,1))*2;       %自由度数，梁单元每个节点有2个自由度，x,y方向位移
f=ones(dof,1)*1e4;             %结构整体外载荷矩阵，整体坐标系下
f_loc=zeros(8,1);              %单元外载荷矩阵，局部坐标系下
u=ones(dof,1)*1e2;             %位移矩阵
K=zeros(dof);                  %总体刚度矩阵

for i=1:n_ele
    k_ele=RectangleElementStiffness(E,miu,t,node(ele(i,2:5),2:3),1);
    K=assemRectangle(K,k_ele,ele(i,2:5));
end

%力边界条件
f(3)=0;          %2节点横向力
f(4)=0;          %2节点垂向力
f(5)=0;      %3节点横向力
f(6)=0;          %3节点垂向力
f(9)=0;      f(10)=0;
f(11)=0;     f(12)=0;
f(15)=0;	f(16)=-100;
f(17)=0;	f(18)=-50;

%位移边界条件
u(1)=0;u(2)=0;
u(7)=0;u(8)=0;
u(13)=0;u(14)=0;

%求解
index=[];      %所需求解自由度的索引
p=[];          %所需求解自由度对应的节点力矩阵
for i=1:dof
    if u(i)~=0
        index=[index,i];
        p=[p;f(i)];
    end
end
u(index)=K(index,index)\p;    %高斯消去
f=K*u;

%单元应力
num_ele=size(ele,1);        %单元数
stress=zeros(num_ele,3);    %应力存储矩阵
x1=node(:,2)+u(1:2:18);     %变形后坐标
y1=node(:,3)+u(2:2:18); 	%这个也记得要改


for i=1:n_ele
    u1=[u(2*ele(i,2)-1);u(2*ele(i,2));u(2*ele(i,3)-1);u(2*ele(i,3));u(2*ele(i,4)-1);u(2*ele(i,4));u(2*ele(i,5)-1);u(2*ele(i,5))];
    stress(i,:)=RectangleElementStress(E,miu,node(ele(i,2:5),2:3),u1,1)';   %单元应力计算
end



%计算单元刚度矩阵
function k_ele=RectangleElementStiffness(E,miu,h,node_ele,p)
%TriangleElementStiffness This function returns the element
% stiffness matrix for a grid
% element with modulus of elasticity E,
% Poission's ratio miu, constant thickness h,
% node_ele the node coordinate of element，plane stress or plane strain option.
% The size of the element stiffness
% matrix is 8 x 8.
syms t s;       
%---------node coordinate------
x1=node_ele(1,1);                y1=node_ele(1,2);
x2=node_ele(2,1);                y2=node_ele(2,2);
x3=node_ele(3,1);                y3=node_ele(3,2);
x4=node_ele(4,1);                y4=node_ele(4,2);
%-------------------------------
%-------shape fuction------------
N1=((1-s)*(1-t))/4;N2=((1+s)*(1-t))/4;
N3=((1+s)*(1+t))/4;N4=((1-s)*(1+t))/4;
%--------------------------------------
%------- Jacobian determinant----------
xc=[x1 x2 x3 x4];
yc=[y1 y2 y3 y4];
J_m=[0 1-t t-s s-1;
     t-1 0 s+1 -s-t;
     s-t -s-1 0 t+1;
     1-s s+t -t-1 0];
J=xc*J_m*yc'/8;
%------------------------------------
%----------gradient function---------
a=(y1*(s-1)+y2*(-1-s)+y3*(1+s)+y4*(1-s))/4;
b=(y1*(t-1)+y2*(1-t)+y3*(1+t)+y4*(-1-t))/4;
c=(x1*(t-1)+x2*(1-t)+x3*(1+t)+x4*(-1-t))/4;
d=(x1*(s-1)+x2*(-1-s)+x3*(1+s)+x4*(1-s))/4;
N1s=diff(N1,s);N1t=diff(N1,t);
N2s=diff(N2,s);N2t=diff(N2,t);
N3s=diff(N3,s);N3t=diff(N3,t);
N4s=diff(N4,s);N4t=diff(N4,t);
B1=[a*N1s-b*N1t 0;
    0 c*N1t-d*N1s;
    c*N1t-d*N1s a*N1s-b*N1t];
B2=[a*N2s-b*N2t 0;
    0 c*N2t-d*N2s;
    c*N2t-d*N2s a*N2s-b*N2t];
B3=[a*N3s-b*N3t 0;
    0 c*N3t-d*N3s;
    c*N3t-d*N3s a*N3s-b*N3t];
B4=[a*N4s-b*N4t 0;
    0 c*N4t-d*N4s;
    c*N4t-d*N4s a*N4s-b*N4t];
B=[B1 B2 B3 B4]/J;
%-------------------------------------
%------strss/strain matrix-------------
if p==1
    D=E/(1-miu^2)*[1 miu 0;
                   miu 1 0;
                   0 0 (1-miu)/2];           %strss/strain matrix for plane stress
elseif p==2
    D=E/(1+miu)/(1-2*miu)*[1-miu miu 0;
                           miu 1-miu 0;
                       0 0 (1-2*miu)/2];     %strss/strain matrix for plane strain
end
BD=B'*D*B*J;
r=int(int(BD,t,-1,1),s,-1,1);
k_ele=double(h*r);    %的单元刚度矩阵
end

%刚度组装函数
function k_t=assemRectangle(k_t,k_ele,node)
%assemRectangle This function assembles the element stiffness
% matrix k of the plane Rectangle element into the global stiffness matrix K.
% This function returns the global stiffness
% matrix K after the element stiffness matrix
% k is assembled.
d(1:2)=2*node(1)-1:2*node(1);d(3:4)=2*node(2)-1:2*node(2);
d(5:6)=2*node(3)-1:2*node(3);d(7:8)=2*node(4)-1:2*node(4);
for ii=1:8
    for jj=1:8
        k_t(d(ii),d(jj))=k_t(d(ii),d(jj))+k_ele(ii,jj);
    end
end
end

%单元应力计算函数
function w=RectangleElementStress(E,miu,node_ele,u1,p)
%RectangleElementStress This function returns the element
% stress matrix for a Rectangle
% element with modulus of elasticity E,
% Poission's ratio miu,
% node_ele the node coordinate of element，plane stress or plane strain option.

syms t s;       %定义自然坐标
%---------node coordinate------
x1=node_ele(1,1);                y1=node_ele(1,2);
x2=node_ele(2,1);                y2=node_ele(2,2);
x3=node_ele(3,1);                y3=node_ele(3,2);
x4=node_ele(4,1);                y4=node_ele(4,2);
%-------------------------------
%-------shape fuction------------
N1=((1-s)*(1-t))/4;N2=((1+s)*(1-t))/4;
N3=((1+s)*(1+t))/4;N4=((1-s)*(1+t))/4;
%--------------------------------------
%------- Jacobian determinant----------
xc=[x1 x2 x3 x4];
yc=[y1 y2 y3 y4];
J_m=[0 1-t t-s s-1;
     t-1 0 s+1 -s-t;
     s-t -s-1 0 t+1;
     1-s s+t -t-1 0];
J=xc*J_m*yc'/8;
%------------------------------------
%----------gradient function---------
a=(y1*(s-1)+y2*(-1-s)+y3*(1+s)+y4*(1-s))/4;
b=(y1*(t-1)+y2*(1-t)+y3*(1+t)+y4*(-1-t))/4;
c=(x1*(t-1)+x2*(1-t)+x3*(1+t)+x4*(-1-t))/4;
d=(x1*(s-1)+x2*(-1-s)+x3*(1+s)+x4*(1-s))/4;
N1s=diff(N1,s);N1t=diff(N1,t);
N2s=diff(N2,s);N2t=diff(N2,t);
N3s=diff(N3,s);N3t=diff(N3,t);
N4s=diff(N4,s);N4t=diff(N4,t);
B1=[a*N1s-b*N1t 0;
    0 c*N1t-d*N1s;
    c*N1t-d*N1s a*N1s-b*N1t];
B2=[a*N2s-b*N2t 0;
    0 c*N2t-d*N2s;
    c*N2t-d*N2s a*N2s-b*N2t];
B3=[a*N3s-b*N3t 0;
    0 c*N3t-d*N3s;
    c*N3t-d*N3s a*N3s-b*N3t];
B4=[a*N4s-b*N4t 0;
    0 c*N4t-d*N4s;
    c*N4t-d*N4s a*N4s-b*N4t];
B=[B1 B2 B3 B4]/J;
%-------------------------------------
%------strss/strain matrix-------------
if p==1
    D=E/(1-miu^2)*[1 miu 0;
                   miu 1 0;
                   0 0 (1-miu)/2];           %strss/strain matrix for plane stress
elseif p==2
    D=E/(1+miu)/(1-2*miu)*[1-miu miu 0;
                           miu 1-miu 0;
                       0 0 (1-2*miu)/2];     %strss/strain matrix for plane strain
end
str=D*B*u1;    %单元内的应力场

wcent = subs(str, [s,t], [0,0]);   %单元中心应力值
w = double(wcent);
end
