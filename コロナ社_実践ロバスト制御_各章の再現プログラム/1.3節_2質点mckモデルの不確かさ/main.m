%% 不確かさをはらむモデル生成，制御を施した場合の応答など確認(PP.18-19)
clear
close all
%% グラフのフォーマット
% fontsize=10.5;
% % GUIのフォント
% set(0, 'defaultUicontrolFontName', 'メイリオ');
% % 軸のフォント
% set(groot, 'defaultAxesFontName','メイリオ');
% % タイトル、注釈などのフォント
% set(groot, 'defaultTextFontName','メイリオ');
% % GUIのフォントサイズ
% set(groot, 'defaultUicontrolFontSize', fontsize);
% % 軸のフォントサイズ
% set(groot, 'defaultAxesFontSize', fontsize);
% % タイトル、注釈などのフォントサイズ
% set(groot, 'defaultTextFontSize', fontsize);
% レジェンドの位置
set(groot, 'defaultLegendLocation', 'northoutside')
% レジェンドの縦横
set(groot, 'defaultLegendOrientation', 'horizontal')
% グラフの太さ
set(groot, 'defaultLineLineWidth', 2)
% デフォルトのグラフの色（背景）
set(groot,'defaultFigureColor','w')
%% 実行1.1
%% ノミナルモデルの定義
s=tf('s'); % ラプラス演算子sの定義
M=1; % 質量
Pn=1/(M*s^2); % ノミナルモデル
%% PD制御器
% 近似微分器を使用
alpha=-0.5;
K1=alpha^2*M-2*alpha*M*s/(0.01*s+1); % コントローラ1(1.9)
alpha=-2.5;
K2=alpha^2*M-2*alpha*M*s/(0.01*s+1); % コントローラ2(1.9)
%% ノミナルモデルに対する応答計算
Tn1=feedback(Pn*K1,1); % Tn1=Pn*K/(1+Pn*K)->相補感度関数，指令値から制御量までの伝達関数に一致
Tn2=feedback(Pn*K2,1); % Tn1=Pn*K/(1+Pn*K)->相補感度関数，指令値から制御量までの伝達関数に一致
%% 図を保存
mkdir('図')
figure(1)
step(Tn1,'--',Tn2,10);
ylim([0 1.4]);
grid
legend('\alpha=-0.5','\alpha=-2.5');
saveas(gcf,'図/ノミナルモデルのステップ応答.fig');
figure(2)
bodemag(Tn1,'--',Tn2);
grid
legend('\alpha=-0.5','\alpha=-2.5');
saveas(gcf,'図/ノミナルモデルのボード線図.fig');

%% 実行1.2
%% 摂動モデルの定義(Robust Control Toolboxが必要)
%% urealを用いた変動パラメータの定義
m1=ureal('m1',0.8,'percent',10); % m1は0.8 kgで±10 %の摂動があると仮定
m2=M-m1; % m1+m2=Mより
k=ureal('k',300,'percent',10); % kは300 N/mで±10 %の摂動があると仮定
c=ureal('c',1,'percent',10); % cは1 Ns/mで±10 %の摂動があると仮定
%% Two mass model集合の定義
P=(c*s+k)/(s^2*(m1*m2*s^2+(m1+m2)*c*s+(m1+m2)*k)); % 摂動込みの制御対象(1.10)
%% モデル集合から50通りのモデルを選択
P=usample(P,50);
%% 図の保存
figure(3)
bodemag(Pn,P,'--',{1e1,1e3}) % 摂動モデルのボード線図
legend('nominal','perturbation');
grid
saveas(gcf,'図/摂動モデルのボード線図.fig')

%% 実行1.3
%% 摂動モデルに対する制御応答計算
T1=feedback(P*K1,1); % T=Pn*K/(1+Pn*K)->相補感度関数，指令値から制御量までの伝達関数に一致
T2=feedback(P*K2,1); % T=Pn*K/(1+Pn*K)->相補感度関数，指令値から制御量までの伝達関数に一致
figure(4)
step(T1,'--',T2,10); % ステップ応答
ylim([0 1.4]);
grid
legend('\alpha=-0.5','\alpha=-2.5');
saveas(gcf,'図/摂動モデルのステップ応答.fig');
figure(5)
nyquist(P*K1) % ナイキスト線図(当然開ループ伝達関数が引数です)
axis([-1.5 0.5 -1 1]);
legend('\alpha=-0.5')
saveas(gcf,'図/摂動モデルのナイキスト線図_alpha=-0.5.fig');
figure(6)
nyquist(P*K2) % ナイキスト線図(当然開ループ伝達関数が引数です)
axis([-1.5 0.5 -1 1]);
legend('\alpha=-2.5')
saveas(gcf,'図/摂動モデルのナイキスト線図_alpha=-2.5.fig');

%% 結論
% 不確かさがある場合，ゲインをでかくすればノミナル性能は向上するが，ロバスト安定性は損なわれていくのだ
% ロバスト制御をmatlabで組みたかったらRobust Control ToolboxとControl System Toolboxが必要