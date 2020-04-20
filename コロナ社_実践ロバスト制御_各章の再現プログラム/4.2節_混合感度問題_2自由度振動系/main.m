%% 混合感度問題(PP.53-69)
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
%% 制御対象の定義
%% パラメータ決定
m1=0.8;
m2=0.2;
k1=100;
k2=ureal('k2',300,'percent',20);
c1=1;
c2=0.3;
Ks=100; % 力定数(N/V)
%% 運動方程式のM,C,Kの定義
M=[m1 0;0 m2];
C=[c1+c2 -c2;-c2 c2];
K=[k1+k2 -k2;-k2 k2];
F=[Ks;0];
%% 状態空間表現へ（1入力1出力になっています）
Ap=[zeros(2,2) eye(2);-M\K -M\C];
Bp=[zeros(2,1);M\F];
Cp=[0 1 0 0]; % 質点2の変位だけが計測可能
Dp=0;
%% 制御対象の定義
P=ss(Ap,Bp,Cp,Dp);
mkdir('図');
figure(1)
bode(P,{1e0,1e2}) % ボード線図
grid
legend('摂動を含む系のボード線図（電圧V->変位x2）')
saveas(gcf,'図/摂動を含む系のボード線図')
% 摂動によって2次共振周波数にぶれ
%% 乗法的摂動としての見積もり及び重みWmの決定
w=logspace(0,3,100); % 周波数ベクトルの定義
P_g=ufrd(P,w); % 周波数応答
Dm_g=(P_g-P_g.nominal)/P_g.nominal; % 乗法的摂動の計算
% 以下でWmを決定
s=tf('s');
Wm=3*s^2/(s^2+18*s+45^2); % 乗法的摂動に関する重み
% ゲイン線図のプロット
figure(2)
bodemag(Dm_g,'--',Wm,'r-',w)
grid
legend('\Delta_m','Wm')
% ここで常にWmがDm_gよりも大きければおっけー
saveas(gcf,'図/乗法的摂動と重み関数のゲイン線図.fig');
%% 感度関数に対する重みを考え，混合感度問題へ帰着
%% 各種重み関数の定義
s=tf('s');
Ws=15/(s+1.5e-2); % 感度関数に関する重み
Wt=Wm; % 乗法的摂動に関する重み
Weps=5e-4; % w->yの伝達関数がフルランクになるための工夫
figure(3)
bodemag(Ws,Wt,'--',w)
legend('Ws','Wt');
grid
saveas(gcf,'図/重み関数のゲイン線図');
%% 一般化プラントの導出
Pn=P.nominal; % ノミナルプラント
systemnames='Pn Ws Wt Weps';
inputvar='[w1;w2;u]';
outputvar='[Ws; Wt; Pn+Weps]';
input_to_Pn='[w1-u]';
input_to_Wt='[u]';
input_to_Ws='[w1-u]';
input_to_Weps='[w2]';
G=sysic;
%% 一般化プラントに対してH∞制御系設計
[K,clp,gamma_min,hinf_info]=hinfsyn(G,1,1,'display','on'); % hinfsynでH∞問題を解いていまーす
%% H∞制御器と制御対象のゲイン線図
figure(4)
w=logspace(0,2,100);
bodemag(K,Pn,'--',w);
legend('K','P.nominal')
grid
saveas(gcf,'図/H∞制御器と制御対象のゲイン線図');
%% 伝達関数が仕様を満たしていることを確認など
% ここで実線部が点線部を下回ることを確認(だめならWtとかWsの設計をやり直しでございます)
T=P*K/(1+P*K); % 相補感度関数
S=1/(1+P*K); % 感度関数
M=P/(1+P*K); % 外乱dから出力yまでの伝達関数
figure(5)
bodemag(T.nominal,'r',S.nominal,'b',1/Wt,'r--',1/Ws,'b--');
legend('T','S','1/Wt','1/Ws');
xlim([10^0 10^2]);
ylim([-30 30]);
grid
saveas(gcf,'図/H∞制御器が設計仕様を満たしているかの確認')
figure(6)
step(T,2);
grid
saveas(gcf,'図/摂動を含むシステムに対するステップ応答')
figure(7)
subplot(211),impulse(M,2)% 外乱d->yへのインパルス応答
ylim([-10 15])
grid
subplot(212),impulse(Pn,2)% ノミナルモデルのインパルス応答
ylim([-10 15])
grid
saveas(gcf,'図/インパルス外乱が出力に与える影響')
%% 結論
% fig6より，全ての摂動に関しては安定化可能
% fig7より，外乱の影響がろくによくなっていないことがわかる，
% これはP,Kの間で極零相殺を起こすことで性能を上げているから
% 以下はそのことの確認，確かに虚軸に近い極について極零相殺が起こっている
pole(Pn)
zero(K)
% 安定な極零相殺は内部安定性に影響を与えないから設計時に問題にならないんだなあ
% しかし，伝達関数の形を見てわかる通り，外乱d->yの伝達関数(102行目)では極零相殺が起こらずプラントの極がそのまま残る
% このままでは摂動に対するロバスト安定化は可能だが，外乱抑制はできないっていうこと
% じゃあどうするの？→修正混合感度問題に問題をすり替える（しなきゃいけないことは割と変わらないけどイメージはわきにくくなる）