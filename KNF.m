function [varargout]=KNF(varargin)
% KNF  Приводит систему управления к нормальной форме Коши.
% 
%   [A,B]=KNF(n,d) Возвращает матричную форму системы, для системы управления, в
%    которой n,d - полиномы числителя и знаменателя системы управления соответственно.
% 
%   [A,B]=KNF(n,d,h) h - шаг интегрирования
% 
%   [A,B]=KNF(n,d,h,t) Для построения переходной характеристики системы управления 
%   относительно полученных матриц, укажите дополнительный параметр t,
%   задающий время моделирования системы.
% 
%   [A,B]=KNF(n,d,h,t,'euler') В качестве метода численного интегрирования по умолчанию используется 
%   метод Адамса 5-го порядка. Также есть возможность использовать метод
%   Эйлера и Рунге-Кутта, для которых требуется указать дополнительный
%   параметр при вызове функции. 'euler' 'runge'
%   
%   В случае, если система управления содержит неизвестные параметры, перед
%   началом построения вам потребуется ввести их в выплывающем окне.

% Получение числа входных и выходных переменных

ni = nargin;
no = nargout;

if ni<2
 warning('Недостаточно входных аргументов');
 help KNF
 return
elseif ni>6
   warning('Слишком много входных аргументов');
 help KNF
 return 
end

%  =====================================================

n=varargin{1};
d=varargin{2};


% =====================================================================
%Делаем векотра числителя и знаменателя соразмерными
numel_n=numel(n);
numel_d=numel(d);

raz=max(numel_d,numel_n)-min(numel_d,numel_n);
if raz~=0
    if numel_d>numel_n
        n=[zeros(1,raz) n];
        n=n/(d(1));
        d=d/(d(1));
    else
       
        dd=[zeros(1,raz) d]
        d=d/(d(1));
        n=n/(d(1));
    end
end

% ======================= A ===================
if numel(d)==2
    A=-d(end);
else
res1=eye(numel(d)-2);

for i=1:(numel(d)-2)
    A(i,:)=[0 res1(i,:)];
end

kk=numel(d);

for j=1:(numel(d)-1)
Acc(j)=-d(kk-j+1);
end
A=[A;Acc];
end
% ======================= B ===================

f(1)=n(1);
for i=2:(numel(n))
    buf=0;
    for k=1:(i-1)
        n(i-k);
        buf=buf+d(i-k+1)*f(k);
    end
f(i)=n(i)-buf;
end
for j=1:numel(f)-1
   B(j)=f(j+1);
end

% ================= plot ======================
if numel(varargin)>=3
    if eq(isempty(symvar(sym(A))),0) || eq(isempty(symvar(sym(B))),0)
        
        u_s=unique([symvar(sym(A)) symvar(sym(B))]);
        pro=cell(numel(u_s),1);

        if numel(varargin)==5
                        for j=1:numel(u_s)
                           pro{j}=char(u_s(j)); 
                        end
                        answer=inputdlg(pro);
                    
        elseif numel(varargin)==6
                        kek=varargin{6};
                        for i=1:numel(kek)
                        answer{i,1}=num2str(kek(i));
                        end

        end

      if sum(isnan(str2double(answer)))~=0
         warning('Введены некорректные значения'); 
         return
      else
%           subs

    A=subs(A,u_s,str2double(answer)');
    B=subs(B,u_s,str2double(answer)');

      end
        

end
    
   h=varargin{3};

    t=0:h:varargin{4};
    
    X(:,1)=zeros(numel(A(:,1)),1);
 
    
  if ni==4
    varargin{5}='adams';
  end

if (strcmp(varargin{5},'dnplot') || strcmp(varargin{5},'euler') || strcmp(varargin{5},'adams') || strcmp(varargin{5},'runge')) && numel(varargin)>=5
    switch varargin{5}
        case 'euler'
            for i=1:numel(t)-1
    X(:,i+1)=X(:,i)+(h)*(A*X(:,i)+B');
            end
            varargout{3}=X;
%             =====================================================
        case 'runge'
                  for i=1:2
    X(:,i+1)=X(:,i)+(h)*(A*X(:,i)+B');
                  end 
      for i=1:numel(t)-3  
    K1=h*(A*X(:,i)+B');
    K2=h*(A*(X(:,i+1)+0.5*K1)+B');
    K3=h*(A*(X(:,i+2)-K1+2*K2)+B');
    X(:,i+3)=(1/6)*(K1+4*K2+K3);
      end
      varargout{3}=X;
%       =============================================================
        case 'adams'  
             for i=1:4
    X(:,i+1)=X(:,i)+(h)*(A*X(:,i)+B');
            end           
 for i=4:numel(t)-1           
X(:,i+1)=X(:,i)+(h/24)*(55*(A*X(:,i)+B')-59*(A*X(:,i-1)+B')+37*(A*X(:,i-2)+B')-9*(A*X(:,i-3)+B'));
 end
varargout{3}=X;
% ========================================================================
case 'dnplot'  
  
             for i=1:4
    X(:,i+1)=X(:,i)+(h)*(A*X(:,i)+B');
            end           
 for i=4:numel(t)-1           
X(:,i+1)=X(:,i)+(h/24)*(55*(A*X(:,i)+B')-59*(A*X(:,i-1)+B')+37*(A*X(:,i-2)+B')-9*(A*X(:,i-3)+B'));
 end
varargout{3}=X;
% ==================================================================
    end
    
elseif numel(varargin)==4
                for i=1:numel(t)-1
    X(:,i+1)=X(:,i)+(h)*(A*X(:,i)+B');
                end
elseif numel(varargin)==5
    warning('Ошибка при выборе метода. Используйте: euler runge adams.');
    return
end

    

    if strcmp(varargin{5},'dnplot')
    else
    figure(1)
    subplot(2,1,1);
    Xw=X(1,:);
    plot(t,Xw);grid on;
    subplot(2,1,2);
    
    
    if exist('u_s')==1
        q1=double(subs(varargin{1},u_s,str2double(answer)'));
        q2=double(subs(varargin{2},u_s,str2double(answer)'));


    step(tf([q1],[q2]),varargout{4});grid on;   
    else
     step(tf([varargin{1}],[varargin{2}]),varargin{4});grid on;   
    end
    
    end
     
end
varargout{1}=A(:,:);
varargout{2}=B(:,:);
% varargout{3}=X(:,:);
% ==================     errors =========================

%  Рассмотреть необходимость правки скрипта для случаяБ когда степень 
%  числителя больше степени знаменателя 

% =======================  fixed =====================
% [a,b]=KNF([300],[1 0 300],10)
% вообще не строит
% проблема возникает в связи с наличием нулевого коэффициента среди
% промежуточной переменной
% /ВР/  полностью перейти на полиномиальтные вычисления

% [a,b]=KNF([300],[1 1 300],10)
% неправильное построение
% /BP/ Матрицы Коши кажутся верными. Проблему решило увеличение шага
% интегрирования. Следовательно проблема была вызвана решателем.(Эйлер)
