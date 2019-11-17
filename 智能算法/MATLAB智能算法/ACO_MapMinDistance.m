function main() 

G=[0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0; 
   0 1 1 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0; 
   0 1 1 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0; 
   0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0; 
   0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0; 
   0 1 1 1 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0; 
   0 1 1 1 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0;
   0 1 1 1 0 0 1 1 1 0 0 1 1 1 0 0 0 0 0 0; 
   0 1 1 1 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0; 
   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
   0 0 0 0 0 0 0 1 1 1 0 1 1 1 0 0 0 0 0 0; 
   0 0 0 0 0 0 0 1 1 1 0 1 1 1 0 0 0 0 0 0; 
   0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0; 
   0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0; 
   1 1 1 1 0 0 0 0 0 0 1 1 1 1 0 1 1 1 1 0; 
   1 1 1 1 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 0; 
   0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 1 1 0; 
   0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 1 1 0; 
   0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 0 0 0 0 0; 
   0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0;];

% G = randi([0,1],10);

% G ����ͼΪ01�������Ϊ1��ʾ�ϰ��� 
MM=size(G,1);     
% Tau ��ʼ��Ϣ�ؾ���
Tau=ones(MM*MM,MM*MM);        
Tau=8.*Tau; 
% ����������ָ���ϳ������ٲ���
K=100;             
% ���ϸ���
M=50;  
% ���·������ʼ��
S=1; 
% ���·����Ŀ�ĵ�
E=MM*MM;  

% Alpha ������Ϣ����Ҫ�̶ȵĲ���
Alpha=1;     
% Beta ��������ʽ������Ҫ�̶ȵĲ���
Beta=7;      
% Rho ��Ϣ������ϵ��
Rho=0.3 ;   
% Q ��Ϣ������ǿ��ϵ�� 
Q=1;  

minkl=inf; 
mink=0; 
minl=0; 
D=G2D(G);
% N��ʾ����Ĺ�ģ�����ظ�����
N=size(D,1); 

% С�������صı߳�
a=1;                     
% ��ֹ�������
Ex=a*(mod(E,MM)-0.5);    
if Ex==-0.5 
    Ex=MM-0.5; 
end 
% ��ֹ��������
Ey=a*(MM+0.5-ceil(E/MM));              
% ����ʽ��Ϣ��ȡΪ��Ŀ����ֱ�߾���ĵ���
Eta=zeros(N);

 %��������ʽ��Ϣ����
for i=1:N 
    ix=a*(mod(i,MM)-0.5); 
    if ix==-0.5 
        ix=MM-0.5; 
    end 
    iy=a*(MM+0.5-ceil(i/MM));  
    if i~=E 
        Eta(i)=1/((ix-Ex)^2+(iy-Ey)^2)^0.5; 
    else 
        Eta(i)=100; 
    end 
end 

% ��ϸ���ṹ�洢ÿһ����ÿһֻ���ϵ�����·��
ROUTES=cell(K,M); 
% �þ���洢ÿһ����ÿһֻ���ϵ�����·�߳���
PL=zeros(K,M); 
                    
% ����K��������ʳ���ÿ���ɳ�Mֻ����
for k=1:K 
    for m=1:M 
        % ״̬��ʼ��
        % ��ǰ�ڵ��ʼ��Ϊ��ʼ��
        W=S;   
        % ����·�߳�ʼ��
        Path=S;  
        % ����·�߳��ȳ�ʼ��
        PLkm=0; 
        % ���ɱ��ʼ��
        TABUkm=ones(N); 
        % �Ѿ��ڳ�ʼ���ˣ����Ҫ�ų�
        TABUkm(S)=0;  
        % �ڽӾ����ʼ��
        DD=D;                 
        % ��һ������ǰ���Ľڵ�
        DW=DD(W,:); 
        DW1=find(DW);
        for j=1:length(DW1) 
            if TABUkm(DW1(j))==0 
                DW(DW1(j))=0; 
            end
        end
        LJD=find(DW); 
        % ��ѡ�ڵ�ĸ���
        Len_LJD=length(LJD);
        % ����δ����ʳ�������������ͬ������ʳֹͣ
        while W~=E&&Len_LJD>=1 
            % ת�ֶķ�ѡ����һ����ô��
            PP=zeros(Len_LJD);
            for i=1:Len_LJD 
                PP(i)=(Tau(W,LJD(i))^Alpha)*((Eta(LJD(i)))^Beta); 
            end
            sumpp=sum(PP); 
            % �������ʷֲ�
            PP=PP/sumpp; 
            Pcum(1)=PP(1);

            for i=2:Len_LJD 
                Pcum(i)=Pcum(i-1)+PP(i); 
            end 
            Select=find(Pcum>=rand); 
            to_visit=LJD(Select(1)); 
            % ״̬���ºͼ�¼
            % ·������
            Path=[Path,to_visit];     
            % ·����������
            PLkm=PLkm+DD(W,to_visit);  
            % �����Ƶ���һ���ڵ�
            W=to_visit;                  

            for kk=1:N 
                if TABUkm(kk)==0 
                    DD(W,kk)=0; 
                    DD(kk,W)=0; 
                end 
            end           
            % �ѷ��ʹ��Ľڵ�ӽ��ɱ���ɾ��
            TABUkm(W)=0;                
            DW=DD(W,:); 
            DW1=find(DW); 
            for j=1:length(DW1) 
                if TABUkm(DW1(j))==0 
                    DW(j)=0; 
                end 
            end 

            LJD=find(DW); 
            % ��ѡ�ڵ�ĸ���
            Len_LJD=length(LJD); 
        end
        % ����ÿһ��ÿһֻ���ϵ���ʳ·�ߺ�·�߳���
        ROUTES{k,m}=Path; 
        if Path(end)==E 
            PL(k,m)=PLkm; 
            if PLkm<minkl 
                mink=k;minl=m;minkl=PLkm; 
            end
        else
            PL(k,m)=0;
        end 
    end 
    % ������Ϣ��
    % ��������ʼ��
    Delta_Tau=zeros(N,N); 
    for m=1:M 
        if PL(k,m)  
            ROUT=ROUTES{k,m}; 
            TS=length(ROUT)-1;%����
            PL_km=PL(k,m); 
            for s=1:TS 
                x=ROUT(s); 
                y=ROUT(s+1); 
                Delta_Tau(x,y)=Delta_Tau(x,y)+Q/PL_km; 
                Delta_Tau(y,x)=Delta_Tau(y,x)+Q/PL_km; 
            end 
        end 
    end 
    Tau=(1-Rho).*Tau+Delta_Tau; % ��Ϣ�ػӷ�һ���֣�������һ����
end 

% ��ͼ
plotif=1; % �Ƿ��ͼ�Ŀ��Ʋ���
% ����������
if plotif==1 
    minPL=zeros(K); 
    for i=1:K 
        PLK=PL(i,:); 
        Nonzero=find(PLK); 
        PLKPLK=PLK(Nonzero); 
        minPL(i)=min(PLKPLK); 
    end      
    
    figure(1) 
    plot(minPL); 
    hold on 
    grid on 
    title('�������߱仯����'); 
    xlabel('��������'); 
    % ������ͼ
    ylabel('��С·������'); 

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
    title('�������˶��켣'); 
    xlabel('����x'); 
    ylabel('����y');
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

    plot(Rx,Ry) 
end

% �������������ͼ
plotif2=0;
if plotif2==1 
    figure(3) 
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

    for k=1:K 
        PLK=PL(k,:); 
        minPLK=min(PLK); 
        pos=find(PLK==minPLK); 
        m=pos(1); 
        ROUT=ROUTES{k,m}; 
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
        
        plot(Rx,Ry) 
        hold on 
    end 
end 
    
function D=G2D(G) 

l=size(G,1); 
D=zeros(l*l,l*l); 
for i=1:l 
    for j=1:l 
        if G(i,j)==0 
            for m=1:l 
                for n=1:l 
                    if G(m,n)==0 
                        im=abs(i-m);jn=abs(j-n); 
                        if im+jn==1||(im==1&&jn==1) 
                            D((i-1)*l+j,(m-1)*l+n)=(im+jn)^0.5; 
                        end 
                    end 
                end 
            end 
        end 
    end 
end

