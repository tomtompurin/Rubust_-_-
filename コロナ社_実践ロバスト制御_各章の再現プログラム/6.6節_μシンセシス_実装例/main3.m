%% 設計 （構造的摂動(実数)＋ロバスト性能）
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
smass_param
J3n=J3;
Kan=Ka;
%% パラメータ摂動の定義
scl=1; % 実数摂動を実数摂動としてみるのでスケールはそのまま
J3=ureal('J3',J3n,'percent',50*scl);
Ka=ureal('Ka',Kan,'percent',10*scl);
%% 摂動モデルPpertとノミナルモデルPn
Ppert=defssmodel(J1,J2,J3,Ka,Kb,D1,D2,D3,Da,Db);
Pn=Ppert.nominal;
%% 相補感度関数に関する重み(修正混合感度問題をそのまま使うため)
Wt=makeweight2(0.1,500,10,0.7);
% P.161, 図6.19の点線とおよそ一致していることを確認
w=logspace(-1,4,500);
figure(1)
bodemag(1/Wt,w)
ylim([-80 30])
title('相補感度関数に関する重みが適切なことの確認')
%% 外乱抑圧特性，観測ノイズに関する重み関数
Wps=tf([1/50 1],[1 1e-3])*500;
Weps=1e-4;
%% 一般化プラントの生成
defgp_mu2
%% D-Kイタレーション
% オプションの定義
dkitopt=dksynOptions(...
    'DisplayWhileAutoIter','on',...
    'MixedMU','on',...
    'NumberOfAutoIterations',7,...
    'FrequencyVector',logspace(1,4,200));
[Kmu3,Gclp,mubnd,dkinfo]=dksyn(G,1,1,dkitopt);
%% 設計結果の確認
mkdir('図/設計3')
load('Kmu1.mat')
load('Kmu2.mat')
figure(1)
bode(Kmu1,'--',Kmu2,'-.',Kmu3);
legend('設計1','設計2','設計3')
grid
saveas(gcf,'図/設計3/コントローラのボード線図')
% 以下の図で重みよりも外乱抑圧特性が下回ればおっけ
figure(2)
bodemag(1/Wps,'--',Ppert/(1+Ppert*Kmu3),w)
legend('1/Wps','P/(1+PK)')
grid
saveas(gcf,'図/設計3/ロバスト性能の確認')
figure(4)
bodemag(1/Wt,'--',Ppert*Kmu3/(1+Ppert*Kmu3),w);
legend('1/Wt','T(=PK/(1+PK))')
grid
saveas(gcf,'図/設計3/相補感度関数')
%% モデルマッチングをする
% 制御対象の相対次数を確認（モデルマッチングのMの次数の決定のため）
tf(Pn)
% Pの相対次数が2なのでモデルマッチングのMの相対次数は2以上
%% 規範モデルM
omega_n=45;
zeta=0.7;
alpha=zeta*omega_n;
s=tf('s');
% Mで定めた規範モデルと同じ応答をさせるような2自由度制御
% omega_nの大きさが大きいほど規範モデルの速応性が高いが，大きくしすぎると，モデルの摂動の影響がある．
M=omega_n^2*alpha/((s^2+2*zeta*omega_n*s+omega_n^2)*(s+alpha));
%% 制御対象Pnの零点を取得
[p,z]=pzmap(Pn);
%% ノッチフィルタの零点とPnの零点が一致するようにノッチフィルタを設計
z_omega=abs(z)
z_zeta=-real(z)./abs(z)
Fn1=notch_f(z_omega(1),z_zeta(1),0.7);
Fn2=notch_f(z_omega(3),z_zeta(3),0.7);
Mt=M;
M=Mt*Fn1*Fn2;
% フィードフォワード制御器
Gff=ss(M/tf(Pn));
Gff=minreal(Gff,1e-05); % 極零相殺を起こす
%% 各ブロックの入出力名を定義
Pn.y='y';
Pn.u='u';
Kmu3.y='ufb';
Kmu3.u='e';
Gff.y='uff';
Gff.u='r';
M.y='yr';
M.u='r';
%% 加算点の定義
sum1=sumblk('u=uff+ufb');
sum2=sumblk('e=yr-y');
%% 各種ブロックの接続
Gtdof=connect(Pn,Kmu3,Gff,M,sum1,sum2,'r',{'y','u'});
% 'r',{'y','u'}は接続したブロックの入出力(CSTの関数です)
%% 2自由度制御による設計の結果
figure(4)
subplot(211)
step(Gtdof(1,1),0.5);
ylim([0 1.2])
title('2自由度制御系の出力')
grid
subplot(212)
step(Gtdof(2,1),0.5);
% ylim([-2 2])
title('2自由度制御系の制御力')
grid
saveas(gcf,'図/設計3/ノッチフィルタありの2自由度制御の結果(ノミナル)')
%% ロバスト性能のチェック
load('Parray')
rng('default');
Ppert=Parray;
Ppert.y='y';
Ppert.u='u';
Gtdof=connect(Ppert,Kmu3,Gff,M,sum1,sum2,'r',{'y','u'});
figure(5)
subplot(211)
step(Gtdof(1,1),0.5);
ylim([0 1.2])
title('2自由度制御系の出力')
grid
subplot(212)
step(Gtdof(2,1),0.5);
ylim([-3 3])
title('2自由度制御系の制御力')
grid
saveas(gcf,'図/設計3/ノッチフィルタありの2自由度制御の結果(ロバスト)')