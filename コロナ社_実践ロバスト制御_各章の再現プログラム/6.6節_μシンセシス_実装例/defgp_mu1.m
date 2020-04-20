%% 一般化プラントを与える（ロバスト性能のみを改善していくようなプラント（図6.15））
systemnames='Ppert Wps Weps'; % H∞ではノミナルモデルを使っていたが，μではモデル集合
inputvar='[d1;d2;u]';
outputvar='[Wps;Ppert+Weps]';
input_to_Ppert='[d1-u]';
input_to_Weps='[d2]';
input_to_Wps='[Ppert+Weps]';
G=sysic;