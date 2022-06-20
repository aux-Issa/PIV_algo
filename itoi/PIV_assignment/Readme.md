PIVデータ（）
のデータ概要

80ppm
 BasePath:/Volumes/HDCZ-UT/itoi_PIV/water/water_test.6uvaasgh/80ppm_experiment.6vaxq7zy
・/velocity: 卒論時に取得したデータ
  BaseFileName:velocity_32px.6vayg6cd.000000.dat
  
・/0217fixed_velocity: 卒論時に取得したデータを加工したもの(座標が違うので間違い)
  BaseFileName:previous_water_Flow_images.6wb7297o.000000.dat
  
・/0415velocity: 卒論時に取得したデータを座標を合わして再解析したもの
  BaseFileName: 0415velocity.6ylneeem.000000.dat
・ /velocity_v2: 壁面の座標を(より壁面よりに範囲を広げた)変更したもの
  BaseFileName: velocity_v2.70e4lubb.000000.dat
  
☆ /velocity_v3: v2にRangeValidationを加えたもの
  BaseFileName: velocity_v3.710rtc9n.000000.dat

water
BasePath: /Volumes/HDCZ-UT/itoi_PIV/water/water_test.6uvaasgh/water_experiment.6vari2x3

・/velocity_v3: 壁面の座標を(より壁面よりに範囲を広げた)変更したもの

・/velocity_v4: v3にRangeValidationを加えたもの
BaseFileName: velocity_v4.710rm3vn.000000.dat
Path: /velocity_v4/velocity_v4.710rm3vn.000000.dat

・実験時に確認すること
速度分布チェック→rmsチェック

ファイル実行順序(sessionに保存される値の関係)
water_algo.m → water.m(sessionに水のデータ保存)
solution_algo.m → solution.m (水と80ppmを対数即に載せる)
          OR
comparision_In_elastozone.m(粘弾性底層内のwater,solution, DNSの比較)

