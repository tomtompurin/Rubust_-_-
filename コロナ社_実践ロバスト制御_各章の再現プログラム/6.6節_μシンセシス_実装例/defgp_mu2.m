%% 一般化プラントを与える（修正混合感度問題（図6.15））
systemnames='Ppert Wps Wt Weps'; % H∞ではノミナルモデルを使っていたが，μではモデル集合
inputvar='[d1;d2;u]';
outputvar='[Wps;Wt;Ppert+Weps]';
input_to_Ppert='[d1-u]';
input_to_Weps='[d2]';
input_to_Wps='[Ppert+Weps]';
input_to_Wt='[u]';
G=sysic;