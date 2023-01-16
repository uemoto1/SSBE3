# SSBE3

作業手順ノート
===========

「Code」ポタンの下にあるDownload zipをクリックするとプログラムパッケージ本体「SSBE3-main.zip」が入手できる。

Wisteriaにインストールする場合:

1. 「SSBE3(略).zip」をworkディレクトリにアップロードする

2. 「unzip SSBE3(略).zip」で圧縮ファイルを展開する

3. 「cd SSBE3(略)/」で展開済みディレクトリに入る

4. 富岳コンピュータ用の設定ファイルをコピーする
cp make.inc.arm-fugaku make.inc

5. コンパイラを有効化する
module load odyssey

6. コンパイルを開始する
make
