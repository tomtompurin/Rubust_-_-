%% 設計1 （非構造的摂動＋ロバスト性能）
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
J3=ureal('J3',J3n,'percent',50);
Ka=ureal('Ka',Kan,'percent',10);
%% 摂動モデルPpertとノミナルモデルPn
Ppert=defssmodel(J1,J2,J3,Ka,Kb,D1,D2,D3,Da,Db);
Pn=Ppert.nominal;
%% 周波数応答
w=logspace(-1,4,500);
rng('default'); % 乱数の初期化
Parray=usample(Ppert,100);
Parray_g=frd(Parray,w);
Pn_g=frd(Pn,w);
mkdir('図/設計1')
figure(1);
bode(Parray_g)
grid
legend('Ppert');
saveas(gcf,'図/設計1/摂動を含むモデル集合');
%% 一応Deltaのサイズだけ確認しておく（最小になるようにdefssmodelを工夫してある）
[G,Delta]=lftdata(Ppert);
size(Delta)
%% 乗法的摂動としての見積もり
Dm_g=(Parray_g-Pn_g)/Pn_g;
%% 乗法的摂動を覆う重み（設計変数）
Wt1=make_wt(28,0.5,30,0.1);
Wt2=make_wt(15,0.7,150,0.7);
Wt=Wt1*Wt1*Wt2*0.35;
Wt_g=frd(Wt,w); % これが重み
% 一応ucoverという関数で重みを自動で作ってもくれるけどハンドメイドの方が性能高いらしいです
% 以下の図で重み(--)が摂動(-)を覆っていれば重みの設計完了
figure(2)
bodemag(Dm_g,Wt_g,'--');
legend('Delta','Wt');
grid
saveas(gcf,'図/設計1/摂動を重みが覆っていることの確認')
%% 不確かさを持つ線形時不変オブジェクトの生成
InputUnc=ultidyn('InputUnc',[1 1]); % 1行1列,H∞ノルムが最大1になる動的摂動を定義
Ppert=Pn*(1+InputUnc*Wt); % 制御モデルの定義
save('Parray','Parray');
%% 感度関数の重みWpsと観測ノイズの重みWepsを定義
Wps=tf([1/50 1],[1 1e-3])*30; % 低周波域でハイゲイン
Weps=1e-4; % 小さく
%% 一般化プラントの生成
defgp_mu1
%% D-Kイタレーション
% オプションの定義
dkitopt=dksynOptions(...
    'DisplayWhileAutoIter','on',...
    'NumberOfAutoIterations',5,...
    'FrequencyVector',logspace(1,4,300));
[Kmu1,clp,bnd]=dksyn(G,1,1,dkitopt);
%% 設計結果の確認
figure(3)
bode(Kmu1)
grid
saveas(gcf,'図/設計1/設計1のコントローラのボード線図')
% 以下の図で重みよりも外乱抑圧特性が下回ればおっけ
figure(4)
bodemag(1/Wps,'--',Ppert/(1+Ppert*Kmu1),w)
legend('1/Wps','P/(1+PK)')
grid
saveas(gcf,'図/設計1/設計1のロバスト性能の確認')
figure(5)
bodemag(Ppert*Kmu1/(1+Ppert*Kmu1),w);
grid
legend('相補感度関数')
saveas(gcf,'図/設計1/設計1の相補感度関数')
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
% フィードフォワード制御器
Gff=ss(M/tf(Pn));
%% 各ブロックの入出力名を定義
Pn.y='y';
Pn.u='u';
Kmu1.y='ufb';
Kmu1.u='e';
Gff.y='uff';
Gff.u='r';
M.y='yr';
M.u='r';
%% 加算点の定義
sum1=sumblk('u=uff+ufb');
sum2=sumblk('e=yr-y');
%% 各種ブロックの接続
Gtdof=connect(Pn,Kmu1,Gff,M,sum1,sum2,'r',{'y','u'});
% 'r',{'y','u'}は接続したブロックの入出力(CSTの関数です)
%% 2自由度制御による設計の結果
figure(6)
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
saveas(gcf,'図/設計1/2自由度制御の結果(ノミナル)')
% この結果から，入力に振動が続く
% これはPが急峻な共振特性を持つからである（Pのボード線図からもわかる）
% Pの共振をキャンセルするように規範モデルMを選択する
% つまりノッチフィルタFを規範モデルMにかける
%% 制御対象Pnの零点を取得
[p,z]=pzmap(Pn);
%% ノッチフィルタの零点とPnの零点が一致するようにノッチフィルタを設計
z_omega=abs(z)
z_zeta=-real(z)./abs(z)
Fn1=notch_f(z_omega(1),z_zeta(1),0.7);
Fn2=notch_f(z_omega(3),z_zeta(3),0.7);
Mt=M;
M=Mt*Fn1*Fn2;
Gff_n=ss(M/tf(Pn));
Gff_n=minreal(Gff_n,1e-05); % 極零相殺を起こす
figure(7)
bodemag(Gff,'--',Gff_n,w)
grid
legend('フィルタなし','フィルタあり')
saveas(gcf,'図/設計1/Gffのゲイン線図')
%% ノッチフィルタ込みでの2自由度制御の結果を見る
%% 各ブロックの入出力名を定義（書かれていないものはフィルタの前後で変更なし）
Gff_n.y='uff';
Gff_n.u='r';
M.y='yr';
M.u='r';
%% 各種ブロックの接続
Gtdof=connect(Pn,Kmu1,Gff_n,M,sum1,sum2,'r',{'y','u'});
% 'r',{'y','u'}は接続したブロックの入出力(CSTの関数です)
%% 2自由度制御による設計の結果
figure(8)
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
saveas(gcf,'図/設計1/ノッチフィルタありの2自由度制御の結果(ノミナル)')
%% ロバスト性能のチェック
rng('default');
Ppert=Parray;
Ppert.y='y';
Ppert.u='u';
Gtdof=connect(Ppert,Kmu1,Gff_n,M,sum1,sum2,'r',{'y','u'});
figure(9)
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
saveas(gcf,'図/設計1/ノッチフィルタありの2自由度制御の結果(ロバスト)')