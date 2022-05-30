PIVデータ（/Volumes/HDCZ-UT/itoi_PIV/water/water_test.6uvaasgh/80ppm_experiment.6vaxq7zy）
のデータ概要

・/velocity: 卒論時に取得したデータ
  velocity_32px.6vayg6cd.000000.dat
・/0217fixed_velocity: 卒論時に取得したデータを加工したもの(座標が違うので間違い)
  previous_water_Flow_images.6wb7297o.000000.dat
・/0415velocity: 卒論時に取得したデータを座標を合わして再解析したもの
  0415velocity.6ylneeem.000000.dat

・実験時に確認すること
速度分布チェック→rmsチェック


ファイル実行順序(sessionに保存される値の関係)
water_algo.m → solution_algo.m 
→ water.m → solution.m
          OR
→ comparision_In_elastozone.m(粘弾性底層内のwater,solution, DNSの比較)