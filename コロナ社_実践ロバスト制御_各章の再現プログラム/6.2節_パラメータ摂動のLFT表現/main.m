%% パラメータ摂動のLFT表現(PP.134-139)
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
%% 実行6.1
%% ノミナル値
m0=1;
k0=100;
c=1;
%% 摂動幅の定義 (ノミナル値の±10 %)
Delta_m=m0*0.1;
Delta_k=k0*0.1;
%% 正規化された実数摂動
delta_1=ureal('delta_1',0);
delta_2=ureal('delta_2',0);
Delta=blkdiag(delta_1,delta_2);
%% 摂動込みのシステムGの定義
A=[0 1; -k0/m0 -c/m0];
B=[0 0 0;-1/m0 -1/m0 1/m0];
C=[-k0/m0*Delta_m -c/m0*Delta_m;Delta_k 0;1 0];
D=[-Delta_m/m0 -Delta_m/m0 Delta_m/m0;0 0 0;0 0 0];
G=ss(A,B,C,D);
%% LFTの計算
P=lft(Delta,G);
%% ボード線図
mkdir('図')
figure(1)
bode(P)
saveas(gcf,'図/パラメータ摂動を含むモデル集合のボード線図')

%% 実行6.2(やりたいことは実は6.1と同一)
%% ノミナル
m0=1;
k0=100;
c=1;
m=ureal('m',m0,'percent',10);
k=ureal('k',k0,'percent',10);
%% 状態空間表現
A=[0 1;-k/m -c/m];
B=[0;1/m];
C=[1 0];
D=0;
P=ss(A,B,C,D);
%% lftへ
[G,Delta,BlkStruc,NormUNC]=lftdata(P);
figure(2)
bode(P)
% ここで，6.1と6.2のDeltaのサイズや，Pの内容をチェックすると，どうやら同じではない
size(Delta)
NormUNC{:}
P
% 同じ構造，同じ摂動でもPが1意に定まらない
% 基本的にDeltaのサイズが小さいほど保守性が低いので小さくするような工夫がいる
% 造ろうと思えば6.1みたいに自力で書けるけど面倒

%% 実行6.7
%% ノミナル
m0=1;
k0=100;
c=1;
m=ureal('m',m0,'percent',10);
k=ureal('k',k0,'percent',10);
%% 状態空間表現（ディスクリプタシステム）
E=diag([1 m]);
A=E\[0 1;-k -c];
B=E\[0;1];
C=[1 0];
D=0;
P=ss(A,B,C,D);
%% lftへ
[G,Delta,BlkStruc,NormUNC]=lftdata(P);
figure(3)
bode(P)
size(Delta)
NormUNC{:}
P
% なんでかディスクリプタシステムを使ってシステムを定義してみるとDeltaのサイズが落ちる（なんでだろう）

%% 結論
% H∞制御では制御系設計時にどの定数が摂動しているか考慮していなかった（非構造的摂動として処理）
% μ設計法では摂動の構造を考慮したり，ロバスト性能を保証する設計が可能です．
% H∞はロバスト性能は結構難しい