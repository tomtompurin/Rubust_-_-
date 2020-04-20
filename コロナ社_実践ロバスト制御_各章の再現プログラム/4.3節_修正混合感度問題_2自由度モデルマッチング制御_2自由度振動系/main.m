%% 修正混合感度問題(PP.69-78)
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
% 凡例の位置
set(groot, 'defaultLegendLocation', 'northoutside')
% 凡例の縦横
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
%% 乗法的摂動としての見積もり及び重みWmの決定（設計者の仕事）
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
%% 重みを考え，修正混合感度問題へ帰着（設計者の仕事）
% やっていることはWsSのH∞ノルムを下げる代わりに，WpsPSのH∞ノルムを下げる問題に変更しただけ
%% 各種重み関数の定義
s=tf('s');
Ws=15/(s+1.5e-2); % 感度関数に関する重み
Wps=0.8*Ws; % こいつが設計変数(こいつは設計結果のガンマが1に近づくように設計するとよろしい)
Wt=Wm; % 乗法的摂動に関する重み
Weps=5e-4; % w->yの伝達関数がフルランクになるための工夫
% Wps*Pが，もともとのWsと比べて，低周波域で小さくなっていないことを確認
figure(3)
bodemag(Ws,'r:',Wps,Wps*P.nominal,'r-.',Wt,'--',w);
legend('Ws','Wps','Wps*P','Wt');
ylim([-60 30])
grid
saveas(gcf,'図/重み関数のゲイン線図');
%% 一般化プラントの導出(P.72,図4.12(b)を再現)
Pn=P.nominal; % ノミナルプラント
systemnames='Pn Wps Wt Weps';
inputvar='[w1;w2;u]';
outputvar='[Wps; Wt; Pn+Weps]';
input_to_Pn='[w1-u]';
input_to_Wt='[u]';
input_to_Wps='[Pn+Weps]';
input_to_Weps='[w2]';
G=sysic;
%% 一般化プラントに対してH∞制御系設計
[K,clp,gamma_min,hinf_info]=hinfsyn(G,1,1,'display','on'); % hinfsynでH∞問題を解いていまーす
% ちなみに引数の1,1は入力数，出力数でした，これで多入力他出力いけそう
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
bodemag(T.nominal,'r',S.nominal,'b',1/Wt,'r--',1/(Wps*P.nominal),'b--');
legend('T','S','1/Wt','1/(Wps*P)');
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
% 修正混合感度問題にすることでロバスト安定化も外乱抑制も可能になっている
% これは極零相殺が起こりにくいように一般化プラントを工夫したからである．
% 設計者の負担を考えると修正混合感度問題の方が面倒だが，
% ロバスト安定化と外乱の抑制を同時にしたいのならばこうするのが一つの手段
% 大きく変わっているのは70行目から92行目で，それ以外はあまりやることが変わらない
% H∞標準問題の仮定を満たす一般化プラントさえ作れてしまえば後はeasy

% ただ，目標値追従性能は修正混合感度問題に変更することで下がってしまっている．
% じゃあどうする？->2自由度制御（一番よさげなのはモデルマッチング制御）
% そもそもH∞ノルムと過渡応答特性に直接的な関係はないのでH∞そのままじゃ無理だね

%% 2自由度制御へ
% 制御対象の相対次数を確認（モデルマッチングのMの次数の決定のため）
tf(P.nominal)
% Pの相対次数が3なのでモデルマッチングのMの相対次数は3以上
%% 規範モデルM
omega_n=20;
zeta=0.7;
alpha=zeta*omega_n;
% Mで定めた規範モデルと同じ応答をさせるような2自由度制御
% omega_nの大きさが大きいほど規範モデルの速応性が高いが，大きくしすぎると，モデルの摂動の影響がある．
M=omega_n^2*alpha/((s^2+2*zeta*omega_n*s+omega_n^2)*(s+alpha));
% フィードフォワード制御器
Gff=ss(M/tf(P.nominal));
%% P,M,Gff,Kを,P80の図4.20(c)のようにつなげる（モデルマッチング型）
%% 各ブロックの入出力名を定義
P.y='y';
P.u='u';
K.y='ufb';
K.u='e';
Gff.y='uff';
Gff.u='r';
M.y='yr';
M.u='r';
%% 加算点の定義
sum1=sumblk('u=uff+ufb');
sum2=sumblk('e=yr-y');
%% 各種ブロックの接続
Gtdof=connect(P,K,Gff,M,sum1,sum2,'r',{'y','u'});
% 'r',{'y','u'}は接続したブロックの入出力(CSTの関数です)

%% 2自由度制御による設計の結果
figure(8)
subplot(211)
step(Gtdof(1,1),1);
ylim([0 1.2])
title('2自由度制御系の出力')
grid
subplot(212)
step(Gtdof(2,1),1);
% ylim([-2 2])
title('2自由度制御系の制御力')
grid
saveas(gcf,'図/2自由度制御の結果')

%% 結論
% ロバスト安定化+ノミナル性能 -> 混合感度問題
% ロバスト安定化+ノミナル性能+外乱抑制 -> 修正混合感度問題
% ロバスト安定化+ノミナル性能+外乱抑制+過渡特性の向上 -> 修正混合感度問題+2自由度モデルマッチング制御

% さらに安定余裕を考慮することが可能で，PP.108-118にやり方が載っています