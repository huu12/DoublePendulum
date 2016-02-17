# DoublePendulum
二重振り子のシミュレーション

## コンパイル
- Mac OS Xの場合 XQuartz が必要です．
- -I や -L は環境によって違う可能性があります．
- Mac OS X 10.11.1, GNU Compiler Collection 5.3.0, XQuartz 2.7.8にて動作確認．
```bash
gcc pendulum.c -I /opt/X11/include -L /opt/X11/lib -lX11 -o pendulum
```

## 実行
```bash
./pendulum
```
```bash
./pendulum l1 l2 m1 m2 theta1 theta2
```
- ln：nつ目の振り子の糸の長さ[m]，m1：nつ目の振り子の質量[kg]，thetan；nつ目の振り子の鉛直下向き方向からの角度[rad]（図参照）
- 引数をなにも取らないと l1=1.4, l2=1.2, m1=1.2, m2=0.4, theta1=3.5, theta2=1.0 とします．
![pendulum](http://i.imgur.com/da7gB5V.png)

## ライセンス
MIT License
