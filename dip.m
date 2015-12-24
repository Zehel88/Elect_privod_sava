function varargout = dip(varargin)
% DIP MATLAB code for dip.fig
%      DIP, by itself, creates a new DIP or raises the existing
%      singleton*.
%
%      H = DIP returns the handle to a new DIP or the handle to
%      the existing singleton*.
%
%      DIP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIP.M with the given input arguments.
%
%      DIP('Property','Value',...) creates a new DIP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dip_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dip_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dip

% Last Modified by GUIDE v2.5 24-Dec-2015 14:14:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dip_OpeningFcn, ...
                   'gui_OutputFcn',  @dip_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end


% --- Executes just before dip is made visible.
function dip_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dip (see VARARGIN)

% Choose default command line output for dip
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% ==========================================================
% preload
set(handles.uitable4,'Visible','off');
set(handles.checkbox1,'Value',0);
load('table.mat');
set(handles.uitable1,'DaTa',[table(1,:);table(2,:);table(3,:)]);

% UIWAIT makes dip wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dip_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --- Run model
function pushbutton2_Callback(hObject, eventdata, handles)
clc
% Подгружаем файл настроек
load('Prefs.mat');
% Подгружаем модель симулинк
load_system('diplom_model.slx');
% Шаг интегрирования
h=str2double(get_param('diplom_model','Fixedstep'));
% Время моделирования
stop_t=str2double(get_param('diplom_model','stoptime'));
t=str2double(get_param('diplom_model','Starttime')):h:stop_t;


% Загружаем передаточную
% Числитель
Num=get_param('diplom_model/tf','Numerator');
% Знаменатель
Den=get_param('diplom_model/tf','Denominator');
% Определяем размерность системы
nnn=numel(Den);
% НУ
x0=zeros(nnn,1);
% Cинтез регулятора
syms s Kp Kd Ki T a
Wo=poly2sym(Num,s)/poly2sym(Den,s);
Wr=Ki/s+Kp+Kd*s;
Wsc=collect(simplify(Wo*Wr/(Wo*Wr+1)),s);
[Numsc Densc]=numden(Wsc);

Numsc_c=rot90(rot90(coeffs(Numsc,s)));

Densc_c=rot90(rot90(coeffs(Densc,s)));

[Asc,Bsc]=KNF(Numsc_c,Densc_c);
Ain=inline(Asc);
Bin=inline(Bsc);

% Цикл для двух методов
for met=1:2
%     Переключатель методов
    switch met
%         Метод алгоритмический
        case 1
%  =======================================================================
            set_param('diplom_model/m2','Gain','1');
            set_param('diplom_model/m1','Gain','0');
%             save_system;
 
%   Задаем начальные параметры регулятора
table=get(handles.uitable1,'DaTa');
% Перезаписали
save('table.mat','table');
% Задаем желаемые параметры регулятора в форме НЧ
params=table;
% Реализация случая с нечеткими параметрами в ОУ
if get(handles.checkbox1,'Value')==1
%    Случай с нечетким ОУ=========================================
% Нечеткие коэф ОУ 
table2=get(handles.uitable4,'DaTa');
koefs = table2;
% Определяем размерность матриц параметров
param_number = size(params,2);
koef_number = size(koefs,2);
% Задаем матрицу индексов для неч пар ОУ
ind =   ones( koef_number ,1);
% Буфферная переменная для параметров ОУ
K = zeros( 3^koef_number, koef_number );
% Буфферная переменная для параметров регулятора
P = zeros( 3^param_number, param_number );

isbreak = 0; j = 1;
while ( isbreak ~= 1 )
    % Формируем матрицу коэф ОУ
    for k = 1:koef_number
        K(j,k) = koefs( ind(k) , k);
    end
    % Индексы матрицы
    isset = 0; k = 1;
    while (isset ~= 1)
        if( k > koef_number )
            isset = 1; isbreak = 1; break
        elseif( ind(k) < 3 )
            ind(k) = ind(k) + 1; isset = 1;
        else
            ind(k) = 1; k = k + 1;
        end
    end
    j = j+1;    
end        



ind =   ones( param_number ,1);
% Cоздадим переменные для заполнения
Y = zeros( 3^koef_number, length(t) );
Y_j = zeros( length(x0), length(t) );  
Y_fuz = zeros( 3 , length(t) );       
W = zeros(3^param_number, 1);         

% Эталон
[Aet,Bet,Xet]=KNF([1],[Prefs{1,2}/2.9 1],h,str2double(get_param('diplom_model','stoptime')),'dnplot');
axes(handles.axes1);
plot(t,Xet);grid on;
xlabel('График эталонной функции');

continue_process = 10;
while( continue_process ~= 0 )

isbreak = 0;
j = 1;
while ( isbreak ~= 1 )
    % Формируем новую матрицу параметров
    for k = 1:param_number
        P(j,k) = params( ind(k) , k);
    end
    % Индексация матрицы
    isset = 0; k = 1;
    while (isset ~= 1)
        if( k > param_number )
            isset = 1; isbreak = 1; break
        elseif( ind(k) < 3 )
            ind(k) = ind(k) + 1; isset = 1;
        else
            ind(k) = 1; k = k + 1;
        end
    end
    j = j+1;
end 




for cur_par = 1:3^param_number

                for cur_koef = 1:3^koef_number
% %                     don't forget to fix initial conditions
% %                     b1 = ( k1 * P(cur_par,2) )/T;
% %                     Y_j(:,1) = ( x0 - [ 0 b1]' ); % add initial conditions
% % 
% %                     for k = 1:(length(t)-1)
% %                         Y_j(:,k+1) = Y_j(:,k) + h * F( t(k), Y_j(:,k), P(cur_par,:), K(cur_koef,:) );
% %                     end


            % [a,b,Y_j]=KNF(Numsc_c,Densc_c,h,stop_t,'dnplot',[P(cur_par,:) K(cur_koef,:)]);
Y_j=Ch_Int(Ain(P(cur_par,1),P(cur_par,2),P(cur_par,3),K(cur_koef,1),K(cur_koef,2)),Bin(P(cur_par,1),P(cur_par,2),P(cur_par,3),K(cur_koef,1),K(cur_koef,2)),h,t);     
% Y_j=Ch_Int(subs(Asc,[Kp Kd Ki T a],[P(cur_par,:) K(cur_koef,:)]),subs(Bsc,[Kp Kd Ki T a],[P(cur_par,:) K(cur_koef,:)]),h,t);

                    Y(cur_koef,:) = Y_j(1,:);  % add solution to the field
                    cur_koef
                 end
    %   2.2) form the core of output
    % don't forget to fix initial conditions
%     Y_j(:,1) = ( x0 - [ 0 b1]' );
%     for k = 1:(length(t)-1)
%         Y_j(:,k+1) = Y_j(:,k) + h * F( t(k), Y_j(:,k), P(cur_par,:), koefs(2,:) );
%     end

% [a,b,Y_j]=KNF(Numsc_c,Densc_c,h,stop_t,'dnplot',[P(cur_par,:) koefs(2,:)]);
Y_j=Ch_Int(Ain(P(cur_par,1),P(cur_par,2),P(cur_par,3),koefs(2,1),koefs(2,2)),Bin(P(cur_par,1),P(cur_par,2),P(cur_par,3),koefs(2,1),koefs(2,2)),h,t);     

    Y_fuz(2,:) = Y_j(1,:);
    %   2.3) form the fuzzy output estimation
    Y_fuz(1,:) = max( Y );
    Y_fuz(3,:) = min( Y );
    %   2.4) form the weightpoint of this combination of parameters
    %   2.4.1) form the L2-norm of core, versus ethalon as sqrt( sum( (y-yj)^2 ) )
    L2_norm = sqrt( sum( (Xet(1,:) - Y_fuz(2, :)).^2 ) );
    %   2.4.2) form the C-norm of core, versus ethalon as max( abs( y-yj) )
    C_norm = max( Xet(1,:) - Y_fuz(2,:) );
    %   2.4.3) form the weight, based on L2-norm of fyzzy borders.
    L2_fuz = sqrt( sum( (Y_fuz(1,:) - Y_fuz(3, :)).^2 ) );
    %   2.4.4) make the weight of current  coefficient to be SUM of found norm-based values
    W(cur_par,1) = L2_norm + C_norm + L2_fuz;
    cur_par
end

[tmp, min_index] = min(W);

%   3') for better understanding - show the etalon and selected outputs together
%   find appropriate solution field
for cur_koef = 1:3^koef_number
%     b1 = ( k1 * P(min_index,2) )/T;
%     Y_j(:,1) = ( x0 - [ 0 b1]' ); % add initial conditions
%     for k = 1:(length(t)-1)
%         Y_j(:,k+1) = Y_j(:,k) + h * F( t(k), Y_j(:,k), P(min_index,:), K(cur_koef,:) );
%     end

Y_j=Ch_Int(Ain(P(min_index,1),P(min_index,2),P(min_index,3),K(cur_koef,1),K(cur_koef,2)),Bin(P(min_index,1),P(min_index,2),P(min_index,3),K(cur_koef,1),K(cur_koef,2)),h,t);     


% Y_j=Ch_Int(subs(Asc,[Kp Kd Ki T a],[P(min_index,:) K(cur_koef,:)]),subs(Bsc,[Kp Kd Ki T a],[P(min_index,:) K(cur_koef,:)]),h,t)
    Y(cur_koef,:) = Y_j(1,:);  % add solution to the field
end
%   form the core of output
% Y_j(:,1) = ( x0 - [ 0 b1]' ); % add initial conditions
% for k = 1:(length(t)-1)
%     Y_j(:,k+1) = Y_j(:,k) + h * F( t(k), Y_j(:,k), P(min_index,:), koefs(2,:) );
% end
% Y_j=Ch_Int(subs(Asc,[Kp Kd Ki T a],[P(min_index,:) koefs(2,:)]),subs(Bsc,[Kp Kd Ki T a],[P(min_index,:) koefs(2,:)]),h,t);
Y_j=Ch_Int(Ain(P(min_index,1),P(min_index,2),P(min_index,3),koefs(2,1),koefs(2,2)),Bin(P(min_index,1),P(min_index,2),P(min_index,3),koefs(2,1),koefs(2,2)),h,t);     


Y_fuz(2,:) = Y_j(1,:);
%   form the fuzzy output estimation
Y_fuz(1,:) = max( Y );
Y_fuz(3,:) = min( Y );

axes(handles.axes2);
plot( t, Xet(1,:), t, Y_fuz(1,:), t, Y_fuz(2,:), t, Y_fuz(3,:) );
grid on;
legend('Ethalon','Higher border','Core','Lower border');

P( min_index, : )
W( min_index, : )



%   5) here we form a new parameters matrix.
for j = 1:param_number
    % in case we found core to be optimal
    if P( min_index, j ) == params(2,j)
        params(1,j) = ( params(1,j) + params(2,j) ) / 2;   % lower support
        params(3,j) = ( params(3,j) + params(2,j) ) / 2;   % higher support
    % in case lower support is found to be optimal
    elseif P( min_index, j ) == params(1,j)
        params(3,j) = params(2,j);     % core becomes higher support
        % lower support stais itself
        params(2,j) = (params(3,j) + params(1,j))/2; % support is considered to be symmethrical
    % in case higher support is found to be optimal
    elseif P( min_index, j ) == params(3,j)
        params(1,j) = params(2,j);    % core becomes lower support
        % higher support remains itself
        params(2,j) = (params(3,j) + params(1,j))/2;  % support is considered to be symmethrical
    end
end
params

continue_process = continue_process - 1;
% if continue_process == 0
%     continue_process = input(' Do you wish to continue? (0 = No, 1 = Yes) :  ');
%     if continue_process ~= 0
%         close all; continue_process = 10;
%     else
%         clear all; close all; continue_process = 0;
%     end
% end



end



else
% Случай с четким ОУ==============================================

end









%  =======================================================================
%             Метод тулбоксовский
        case 2
            set_param('diplom_model/m2','Gain','0');
            set_param('diplom_model/m1','Gain','1');
            save_system('diplom_model.slx');
            

            
    end % end switch met
    
    
    
    
end  % end for met



function editKp_Callback(hObject, eventdata, handles)
% hObject    handle to editKp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editKp as text
%        str2double(get(hObject,'String')) returns contents of editKp as a double


% --- Executes during object creation, after setting all properties.
function editKp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editKp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editKd_Callback(hObject, eventdata, handles)
% hObject    handle to editKd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editKd as text
%        str2double(get(hObject,'String')) returns contents of editKd as a double


% --- Executes during object creation, after setting all properties.
function editKd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editKd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editKi_Callback(hObject, eventdata, handles)
% hObject    handle to editKi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editKi as text
%        str2double(get(hObject,'String')) returns contents of editKi as a double


% --- Executes during object creation, after setting all properties.
function editKi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editKi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)

if get(handles.checkbox1,'Value')==1
    % on opening
prompt = {'Введите кол-во неч. параметров ОУ:'};
dlg_title = 'Input';
def = {'2'};
res = str2double(inputdlg(prompt,dlg_title,1,def));


T=ones(3,res);
set(handles.uitable4,'DaTa',T);
set(handles.uitable4,'Visible','On');

else
    set(handles.uitable4,'Visible','Off');
end




% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
  set_param('diplom_model/m2','Gain','0');
            set_param('diplom_model/m1','Gain','1');
            load_system('diplom_model.slx');
          set_param('diplom_model/tf','Denominator','[1 1 1]');
          sim('diplom_model.slx')
          
          fres=load('fr_d_m.mat');
          
%           b=fres();
          fres
          
          axes(handles.axes3);
          plot(fres.fr_d_m(1,:),fres.fr_d_m(2,:));

          
              
            


% ----------------Запуск модели симулинк (форма меню)-------------
function Untitled_1_Callback(hObject, eventdata, handles)
%% Открытие модели симелинк 
open('diplom_model.slx');





% -----------------Запуск моделирования (форма меню)---------------
function Untitled_2_Callback(hObject, eventdata, handles)








% -----------Настройки (форма меню)------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
load('Prefs');
prompt = {'Введите время управления:'};
dlg_title = 'Input';
def = {'2'};
res = str2double(inputdlg(prompt,dlg_title,1,def));
Prefs{1,2}=res;
save('Prefs.mat','Prefs');
