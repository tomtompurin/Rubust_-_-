%% 乗法的摂動・加法的摂動の記述(PP.40-46)
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

%% 実行3.1 乗法的摂動
s=tf('s');
Pn=1/(s+1); % ノミナル
Wm=2*s/(s+10); % 乗法的摂動のゲイン特性(当然高周波数ほどでけえ)
delta=ultidyn('delta',[1 1],'SampleStateDim',4);
% ultidyn()->動的摂動を生成するコマンド
% 'delta'->摂動の名前
% [1 1]->摂動のサイズ(1行1列)，伝達関数の入出力のサイズに合わせるはず
% 'SampleStateDim'->摂動の次数
P=(1+Wm*delta)*Pn; % 動的摂動を含むモデル集合の作成
P=usample(P,50); % モデル集合から50個引っ張ってくる
w=logspace(-2,2,100);
mkdir('図')
figure(1)
bodemag(Pn,'--',P,w);
grid
legend('nominal','perturbation');
saveas(gcf,'図/乗法的摂動もつモデル集合のゲイン線図')

%% 実行3.2 あるパラメータが20 %摂動したときの加法的摂動と情報的摂動の見積もり（大事）
% 2次遅れ系のωとζが20 %ずれた時の加法的摂動・情報的摂動ってどれくらいなんだろね
omega_n=ureal('omega',1,'percent',20);
zeta=ureal('zeta',0.1,'percent',20);
%伝達関数の定義
s=tf('s');
P=omega_n^2/(s^2+2*zeta*omega_n*s+omega_n^2);
% 周波数応答の計算
w=logspace(-2,2,100);
P_g=ufrd(P,w); % 摂動を含む場合はfrdではなくufrdをつかうんだなあ
% 周波数応答を求めている，実際のモデルでも周波数応答だけ得られていることが多いからそれは既知としたい
% また，今回みたいにどの定数がどのくらい摂動するかでも見積れる

% 乗法的摂動
Dm_g=(P_g-P_g.nominal)/P_g.nominal;
% 加法的摂動
Da_g=P_g-P_g.nominal;
% ゲイン線図のプロット
figure(2)
bodemag(P_g,Dm_g,'--');
legend('P','\Delta_m');
grid
saveas(gcf,'図/乗法的摂動のゲイン線図');
figure(3)
bodemag(P_g,Da_g,'--');
legend('P','\Delta_a');
grid
saveas(gcf,'図/加法的摂動のゲイン線図');

%% まとめ
% 実行3.2のΔを覆うように摂動のゲインの摂動Wを決めるんだなあ
% 加法的摂動か情報的摂動かはどうやって見抜くんだろうなあ？